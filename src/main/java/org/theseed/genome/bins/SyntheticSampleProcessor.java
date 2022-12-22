/**
 *
 */
package org.theseed.genome.bins;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.BinSynthFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.PatricSynthFactory;
import org.theseed.genome.SynthFactory;
import org.theseed.genome.SynthWriter;
import org.theseed.proteins.kmers.reps.RepGenome;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command creates a synthetic sample for hammer testing from a binning directory and/or a PATRIC evaluation
 * result file.  It searches each sample subdirectory for good GTOs and then writes out their contigs.  It then
 * picks random genomes from the evaluation result file that are mostly-good (fail only because of lack of SSU)
 * and writes out their contigs, too.  All genomes chosen are copied to a specified genome output directory
 * for later analysis.
 *
 * The positional parameters are the name of the target representative genome database file (usually "repXXX.ser"),
 * and the name of an output directory to hold the genomes found.  The contigs will be written to the standard output.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	name of the FASTA output file (if not STDOUT)
 *
 * --contigFrac		fraction of contigs that should be included in the output for each genome; the default is 1.0
 * --binDir			name of a binning directory to use for locating binned genomes
 * --evalFile		name of a PATRIC evaluation file from which to harvest mostly-good genomes
 * --max			maximum number of genomes from a single genome source
 * --clear			erase the genome output directory before processing
 *
 * @author Bruce Parrello
 *
 */
public class SyntheticSampleProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SyntheticSampleProcessor.class);
    /** representative-genome database used to generate the hammers */
    private RepGenomeDb repDb;
    /** list of genome sources to process */
    private Collection<SynthFactory> factories;

    // COMMAND-LINE OPTIONS

    /** fraction of contigs to include */
    @Option(name = "--contigFrac", metaVar = "0.90", usage = "fraction of contigs per genome to include")
    private double contigFrac;

    /** output file name (if not STDOUT) */
    @Option(name = "--output", metaVar = "contigs.fa", aliases = { "-o" }, usage = "name of the output FASTA file (if not STDOUT)")
    private File outFile;

    /** TRUE to erase the output directory before beginning */
    @Option(name = "--clear", usage = "if specified, the output genome directory will be erased before processing")
    private boolean clearFlag;

    /** input binning directory */
    @Option(name = "--binDir", metaVar = "binDir",
            usage = "master binning directory, containing sample output in each subdirectory, for harvesting genomes (optional)")
    private File binDir;

    /** input evaluation file */
    @Option(name = "--evalFile", metaVar = "patric.sort.tbl",
            usage = "evaluation output file, containing PATRIC genome evaluation results, for harvesting genomes (optional)")
    private File evalFile;

    /** maximum number of genomes to use from a source */
    @Option(name = "--max", metaVar = "500", usage = "maximum number of genomes to select from a single source")
    private int maxGenomes;

    /** representative-genome database file */
    @Argument(index = 0, metaVar = "repDbFile.ser", usage = "name of the repgen definition file for the repgen set to use for distances",
            required = true)
    private File repDbFile;

    /** output cache directory */
    @Argument(index = 1, metaVar = "outDir", usage = "genome output directory, later usable as a genome source")
    private File outDir;

    @Override
    protected void setDefaults() {
        this.contigFrac = 1.0;
        this.maxGenomes = 1000;
        this.outFile = null;
        this.binDir = null;
        this.evalFile = null;
        this.clearFlag = false;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Insure we have at least one valid source location.
        if (this.evalFile == null && this.binDir == null)
            throw new ParseFailureException("At least one source location specified.");
        // Verify the tuning numbers.
        if (this.contigFrac <= 0.0 || this.contigFrac > 1.0)
            throw new ParseFailureException("Contig fraction must be between 0 and 1.");
        if (this.maxGenomes < 1)
            throw new ParseFailureException("Maximum-genome threshold must be at least 1.");
        // Insure we have an output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating genome output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else if (this.clearFlag) {
            log.info("Erasing genome output directory {}.", this.outDir);
            FileUtils.cleanDirectory(this.outDir);
        } else
            log.info("Selected genomes will be cached in {}.", this.outDir);
        // Load the repgen database.
        if (! this.repDbFile.canRead())
            throw new FileNotFoundException("RepGen database file " + this.repDbFile + " is not found or unreadable.");
        log.info("Loading representative-genome data from {}.", this.repDbFile);
        this.repDb = RepGenomeDb.load(this.repDbFile);
        // We will build the synthetic-genome factories here.
        this.factories = new ArrayList<SynthFactory>();
        if (this.evalFile != null) {
            if (! this.evalFile.canRead())
                throw new FileNotFoundException("Evaluation results file is not found or unreadable.");
            log.info("Creating PATRIC genome factory using {}.", this.evalFile);
            var factory = new PatricSynthFactory(this.evalFile, this.maxGenomes);
            this.factories.add(factory);
        }
        if (this.binDir != null) {
            if (! this.binDir.isDirectory())
                throw new FileNotFoundException("Binning master directory is not found or unreadable.");
            log.info("Creating binning genome factory using {}.", this.binDir);
            var factory = new BinSynthFactory(this.binDir, this.maxGenomes);
            this.factories.add(factory);
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Set up the output stream.
        try (SynthWriter outWriter = new SynthWriter(this.outFile, this.contigFrac)) {
            // Loop through the factories.
            int gCount = 0;
            for (SynthFactory factory : this.factories) {
                // Loop through the factory output.
                for (Genome gto = factory.getNext(); gto != null; gto = factory.getNext()) {
                    // Now we must find the closest representative.
                    RepGenome gtoSeed = this.repDb.getSeedProtein(gto);
                    if (gtoSeed == null)
                        log.info("No seed protein in {}.", gto);
                    else {
                        gCount++;
                        var repInfo = this.repDb.findClosest(gtoSeed);
                        // Get the important data about this representative.
                        String repId = repInfo.getGenomeId();
                        double dist = repInfo.getDistance();
                        // Output the genome.
                        outWriter.writeGenome(gto, repId, dist);
                        File outFile = new File(this.outDir, gto.getId() + ".gto");
                        // Save it to the genome cache.
                        gto.save(outFile);
                    }
                }
            }
            log.info("All done.  {} genomes output.", gCount);
        }
    }

}
