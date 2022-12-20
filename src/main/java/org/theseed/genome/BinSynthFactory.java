/**
 *
 */
package org.theseed.genome;

import java.io.File;
import java.io.FileFilter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This class searches through a binning master directory to find eligible genomes for synthetic samples.
 * We look at the output GTOs in each sample and return the ones that are mostly good.
 *
 * @author Bruce Parrello
 *
 */
public class BinSynthFactory extends SynthFactory {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinSynthFactory.class);
    /** collection of GTO files to scan */
    private List<File> gtoFiles;
    /** iterator through the GTO file list */
    private Iterator<File> gtoIter;
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

    /**
     * Construct a bin-based synthetic sample factory.
     *
     * @param binDir	master binning directory
     * @param max		maximum number of genomes to return
     */
    public BinSynthFactory(File binDir, int max) {
        super(max);
        // Get all the sample directories.
        File[] subDirs = binDir.listFiles(DIR_FILTER);
        log.info("{} subdirectories found in binning directory {}.", subDirs.length, binDir);
        // Create a collection of all the GTO files.
        this.gtoFiles = Arrays.stream(subDirs).flatMap(x -> Arrays.stream(x.listFiles(BIN_GTO_FILTER)))
                .collect(Collectors.toList());
        log.info("{} bin genomes found in {} samples.", gtoFiles.size(), subDirs.length);
        // Shuffle it to randomize our selection.
        Collections.shuffle(this.gtoFiles);
        // Create an iterator to use to run through it.
        this.gtoIter = this.gtoFiles.iterator();
    }

    @Override
    protected Genome getNextGenome() throws IOException {
        // Loop until we find a new genome or run out of them./
        Genome retVal = null;
        while (this.gtoIter.hasNext() && retVal == null) {
            File gtoFile = this.gtoIter.next();
            log.debug("Processing bin {}.", gtoFile);
            // Load the genome for this bin.
            Genome gto = new Genome(gtoFile);
            JsonObject quality = gto.getQuality();
            if (quality.getBooleanOrDefault(QualityKeys.MOSTLY_GOOD)) {
                log.info("Bin {} selected for output.", gtoFile);
                retVal = gto;
            }
        }
        return retVal;
    }

}
