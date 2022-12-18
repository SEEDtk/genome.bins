/**
 *
 */
package org.theseed.genome.bins;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.QualityKeys;
import org.theseed.genome.SynthWriter;
import org.theseed.proteins.kmers.reps.RepGenome;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This command creates a synthetic sample for hammer testing from a binning directory.  It searches each sample
 * subdirectory for good GTOs and then writes out their contigs.
 *
 * The positional parameters are the name of the target representative genome database file (usually "repXXX.ser"),
 * the name of the input binning subdirectory, and the name of an output directory to hold the genomes found.  The
 * contigs will be written to the standard output.
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
public class BinSyntheticSampleProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinSyntheticSampleProcessor.class);
    /** representative-genome database used to generate the hammers */
    private RepGenomeDb repDb;
    /** match pattern for bin output GTOs */
    private static final Pattern BIN_GTO_PATTERN = Pattern.compile("bin\\.\\d+.\\d+\\.gto");
    /** file filter for binning subdirectories */
    private static final FileFilter DIR_FILTER = new FileFilter() {
        @Override
        public boolean accept(File pathname) {
            return pathname.isDirectory();
        }
    };
    /** file filter for bin output genomes */
    private static final FilenameFilter BIN_GTO_FILTER = new FilenameFilter() {
        @Override
        public boolean accept(File dir, String name) {
            return BIN_GTO_PATTERN.matcher(name).matches();
        }
    };

    // COMMAND-LINE OPTIONS

    /** fraction of contigs to include */
    @Option(name = "--contigFrac", metaVar = "0.90", usage = "fraction of contigs per genome to include")
    private double contigFrac;

    /** output file name (if not STDOUT) */
    @Option(name = "--output", metaVar = "contigs.fa", aliases = { "-o" }, usage = "name of the output FASTA file (if not STDOUT)")
    private File outFile;

    /** representative-genome database file */
    @Argument(index = 0, metaVar = "repDbFile.ser", usage = "name of the repgen definition file for the repgen set to use for distances")
    private File repDbFile;

    /** input binning directory */
    @Argument(index = 1, metaVar = "binDir", usage = "binning input directory, containing sample output in each subdirectory")
    private File binDir;

    /** output cache directory */
    @Argument(index = 2, metaVar = "outDir", usage = "genome output directory, later usable as a genome source")
    private File outDir;

    @Override
    protected void setDefaults() {
        this.contigFrac = 1.0;
        this.outFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify the contig fraction.
        if (this.contigFrac <= 0.0 || this.contigFrac > 1.0)
            throw new ParseFailureException("Contig fraction must be between 0 and 1.");
        // Insure we have an output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating genome output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        }
        // Load the repgen database.
        if (! this.repDbFile.canRead())
            throw new FileNotFoundException("RepGen database file " + this.repDbFile + " is not found or unreadable.");
        log.info("Loading representative-genome data from {}.", this.repDbFile);
        this.repDb = RepGenomeDb.load(this.repDbFile);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Get all the sample directories.
        File[] subDirs = this.binDir.listFiles(DIR_FILTER);
        log.info("{} subdirectories found in binning directory {}.", subDirs.length, this.binDir);
        // Create a collection of all the GTO files.
        Collection<File> gtoFiles = Arrays.stream(subDirs).flatMap(x -> Arrays.stream(x.listFiles(BIN_GTO_FILTER)))
                .collect(Collectors.toList());
        log.info("{} bin genomes found in {} samples.", gtoFiles.size(), subDirs.length);
        // Set up the output stream.
        try (SynthWriter outWriter = new SynthWriter(this.outFile, this.contigFrac)) {
            // Loop through the bin files.
            int gCount = 0;
            for (File gtoFile : gtoFiles) {
                gCount++;
                log.info("Processing bin {} of {}: {}.", gCount, gtoFiles.size(), gtoFile);
                // Load the genome for this bin.
                Genome gto = new Genome(gtoFile);
                JsonObject quality = gto.getQuality();
                boolean goodFlag = quality.getBooleanOrDefault(QualityKeys.MOSTLY_GOOD);
                if (! goodFlag)
                    log.info("Skipping bad bin {}.", gto);
                else {
                    // Now we must find the closest representative.
                    RepGenome gtoSeed = this.repDb.getSeedProtein(gto);
                    if (gtoSeed == null)
                        log.info("No seed protein in {}.", gto);
                    else {
                        var repInfo = this.repDb.findClosest(gtoSeed);
                        // Get the important data about this representative.
                        String repId = repInfo.getGenomeId();
                        double dist = repInfo.getDistance();
                        // Output the genome.
                        outWriter.writeGenome(gto, repId, dist);
                        File outFile = new File(this.outDir, gto.getId() + ".gto");
                        gto.save(outFile);
                    }
                }
            }
            log.info("All done.  {} bins examined.", gCount);
        }
    }

}
