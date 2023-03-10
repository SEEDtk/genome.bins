/**
 *
 */
package org.theseed.genome.bins;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.SynthWriter;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command re-purposes a synthetic sample genome directory for a different repgen set.  The positional
 * parameters are the name of the target representative genome database (usually "repXX.ser") and the name
 * of the directory containing the genomes to be repurposed.  For each genome, we find the closest
 * representative and rewrite the genome contigs with corrected information to an output file.
 *
 * The output FASTA file will be written to the standard output.  This may be modified using the "--output"
 * option.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	name of the FASTA output file (if not STDOUT)
 *
 * --contigFrac		fraction of contigs that should be included in the output for each genome; the default is 1.0
 *
 * @author Bruce Parrello
 *
 */
public class SyntheticSampleRewriteProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SyntheticSampleProcessor.class);
    /** representative-genome database used to generate the hammers */
    private RepGenomeDb repDb;
    /** input genome source */
    private GenomeSource genomes;

    // COMMAND-LINE OPTIONS

    /** fraction of contigs to include */
    @Option(name = "--contigFrac", metaVar = "0.90", usage = "fraction of contigs per genome to include")
    private double contigFrac;

    /** output file name (if not STDOUT) */
    @Option(name = "--output", metaVar = "contigs.fa", aliases = { "-o" }, usage = "name of the output FASTA file (if not STDOUT)")
    private File outFile;

    /** representative-genome database file */
    @Argument(index = 0, metaVar = "repDbFile.ser", usage = "name of the repgen definition file for the repgen set to use for distances",
            required = true)
    private File repDbFile;

    /** output cache directory */
    @Argument(index = 1, metaVar = "inDir", usage = "genome input directory")
    private File inDir;

    @Override
    protected void setDefaults() {
        this.contigFrac = 1.0;
        this.outFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Insure the contig fraction is valid.
        if (this.contigFrac <= 0.0 || this.contigFrac > 1.0)
            throw new ParseFailureException("Contig fraction must be between 0 and 1.");
        // Load the repgen database.
        if (! this.repDbFile.canRead())
            throw new FileNotFoundException("Repgen database file " + this.repDbFile + " not found or unreadable.");
        log.info("Loading repgen database from {}.", this.repDbFile);
        this.repDb = RepGenomeDb.load(this.repDbFile);
        // Connect to the input genome source.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input genome directory " + this.inDir + " is not found or invalid.");
        log.info("Connecting to genome input directory {}.", this.inDir);
        this.genomes = GenomeSource.Type.DIR.create(this.inDir);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Set up the output writer.
        try (SynthWriter writer = new SynthWriter(this.outFile, this.contigFrac)) {
            // Loop through the genomes.
            log.info("{} genomes will be processed.", this.genomes.size());
            int count = 0;
            int outCount = 0;
            for (Genome genome : this.genomes) {
                count++;
                log.info("Processing genome {}: {}.", count, genome);
                // Find the representative.
                var seed = this.repDb.getSeedProtein(genome);
                if (seed == null)
                    log.error("ERROR: No seed protein found in {}.", genome);
                else {
                    var rep = this.repDb.findClosest(seed);
                    if (rep.getSimilarity() < this.repDb.getThreshold())
                        log.warn("WARNING: no close neighbor found for {}.", genome);
                    else {
                        writer.writeGenome(genome, rep.getGenomeId(), rep.getDistance());
                        outCount++;
                    }
                }
            }
            log.info("{} genomes written to output.", outCount);
        }
    }

}
