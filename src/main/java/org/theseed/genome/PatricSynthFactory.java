/**
 *
 */
package org.theseed.genome;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Genome;

/**
 * This synthetic-sample factory reads a PATRIC evaluation result file and selects mostly-good genomes that are not
 * good, that is genomes that would be considered good if they had a quality SSU.  The resulting genomes are
 * downloaded from PATRIC.
 *
 * @author Bruce Parrello
 *
 */
public class PatricSynthFactory extends SynthFactory {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(PatricSynthFactory.class);
    /** collection of IDs for genomes to return */
    private Collection<String> genomeIDs;
    /** PATRIC connection for downloading the genomes */
    private P3Connection p3;
    /** iterator through the genome ID collection */
    private Iterator<String> genomeIter;

    /**
     * Construct a PATRIC-genome synthetic-sample factory.
     *
     * @param evalFile	PATRIC genome evaluation file
     * @param max		maximum number of genomes to return
     *
     * @throws IOException
     */
    public PatricSynthFactory(File evalFile, int max) throws IOException {
        super(max);
        // Create the genome ID collection.
        Set<String> eligibleSet = new HashSet<String>((int) (evalFile.length() / 400));
        // We need some counters for tracing.
        int lineCount = 0;
        int skipCount = 0;
        int badCount = 0;
        // We need to read the incoming file to find mostly-good genomes.  Our criterion for this is 90% complete,
        // less than 10% contamination, 80% fine consistency, less than 30% hypotheticals, and a good PheS.  This
        // is a little different from the standard, but it is designed to insure we have useful roles in the genome.
        try (TabbedLineReader inStream = new TabbedLineReader(evalFile)) {
            int idColIdx = inStream.findField("Genome");
            int seedColIdx = inStream.findField("Good Seed");
            int goodColIdx = inStream.findField("Good");
            int hypoColIdx = inStream.findField("Hypothetical");
            int consisColIdx = inStream.findField("Fine");
            int compColIdx = inStream.findField("Completeness");
            int contamColIdx = inStream.findField("Contamination");
            log.info("Scanning {} for mostly-good genomes.", evalFile);
            for (var line : inStream) {
                lineCount++;
                // Skip over genomes that are good or have no seed protein.
                if (line.getFlag(goodColIdx) || ! line.getFlag(seedColIdx))
                    skipCount++;
                else {
                    // Check the numerical constraints.
                    if (line.getDouble(hypoColIdx) > 30.0 || line.getDouble(consisColIdx) < 80.0 || line.getDouble(contamColIdx) > 10.0
                            || line.getDouble(compColIdx) < 90.0)
                        badCount++;
                    else {
                        // Here the genome is good enough.
                        eligibleSet.add(line.get(idColIdx));
                    }
                }
                if (log.isInfoEnabled() && lineCount % 5000 == 0)
                    log.info("{} genomes processed:  {} skipped, {} bad, {} kept.", lineCount, skipCount, badCount, eligibleSet.size());
            }
            log.info("Randomizing selection of {} eligible genomes.", eligibleSet.size());
            this.genomeIDs = this.selectGenomes(eligibleSet);
        }
        log.info("{} genomes processed:  {} skipped, {} bad, {} of {} selected.", lineCount, skipCount, badCount,
                this.genomeIDs.size(), eligibleSet.size());
        // Create the iterator.
        this.genomeIter = this.genomeIDs.iterator();
        // Connect to PATRIC.
        this.p3 = new P3Connection();
    }

    @Override
    protected Genome getNextGenome() throws IOException {
        Genome retVal = null;
        if (this.genomeIter.hasNext()) {
            String genomeID = this.genomeIter.next();
            log.info("Downloading {}.", genomeID);
            retVal = P3Genome.load(this.p3, genomeID, P3Genome.Details.FULL);
        }
        return retVal;
    }

}
