/**
 *
 */
package org.theseed.genome.bins;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.counters.Spacer;
import org.theseed.genome.Contig;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Genome;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseInputProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command creates a synthetic sample based around a repgen set.  The standard input should contain
 * a repgen list file, which contains a set of genomes (ID in "genome_id", name in "genome_name") along
 * with the closest representative (ID in "rep_id", name in "rep_name") and the seed protein distance ("distance").
 * We will select a specified number of genomes, distributed as evenly as possible across the repgen sets,
 * and for each one we will output its contigs.  None of the representative genomes themselves
 * will be selected.
 *
 * The contigs output will have the genome ID and contig ID as the contig label, with an intervening colon (:).
 * This can lead to pretty ugly labels, e.g. "100226.1:100226.1.con.001"; however, the idea is to have something
 * that can be analyzed easily by a computer program regardless of the contig ID format.  The comment of the
 * contig sequence will be the genome name, the repgen ID, and the distance, all separated by tabs.
 *
 * There are no positional parameters.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	name of the input repgen list file (if not STDIN)
 * -o	name of the FASTA output file (if not STDOUT)
 *
 * --min			minimum number of genomes to output (there will be at least one per representative with a neighborhood)
 * --contigFrac		fraction of contigs that should be included in the output for each genome; the default is 1.0
 *
 * @author Bruce Parrello
 *
 */
public class SyntheticSampleProcessor extends BaseInputProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SyntheticSampleProcessor.class);
    /** map of representative genomes to neighbor lists */
    private Map<String, List<Neighbor>> neighborhoodMap;
    /** neighbor genome ID column index */
    private int genomeIdColIdx;
    /** neighbor genome name column index */
    private int genomeNameColIdx;
    /** representative ID column index */
    private int repIdColIdx;
    /** distance column index */
    private int distanceColIdx;
    /** connection to BV-BRC */
    private P3Connection p3;
    /** empty neighbor list to use for genome sorter */
    private static final List<Neighbor> NO_NEIGHBORS = Collections.emptyList();

    // COMMAND-LINE OPTIONS

    /** minimum number of genomes to output */
    @Option(name = "--min", metaVar = "200", usage = "minimum number of genomes to output")
    private int minGenomes;

    /** fraction of contigs to include */
    @Option(name = "--contigFrac", metaVar = "0.90", usage = "fraction of contigs per genome to include")
    private double contigFrac;

    /** output file name (if not STDOUT) */
    @Option(name = "--output", metaVar = "contigs.fa", aliases = { "-o" }, usage = "name of the output FASTA file (if not STDOUT)")
    private File outFile;

    /**
     * This is a simple utility class that contains a neighboring genome, its name, and its distance to the
     * representative.  It is sorted by distance (closest to furthest).
     */
    protected static class Neighbor implements Comparable<Neighbor> {

        /** ID of the neighbor genome */
        private String genomeID;
        /** name of the neighbor genome */
        private String name;
        /** distance to the representative */
        private double distance;

        /**
         * Create a neighbor descriptor.
         *
         * @param genome_id		ID of the neighbor genome
         * @param genome_name	name of the neighbor genome
         * @param dist			distance to the representative
         */
        protected Neighbor(String genome_id, String genome_name, double dist) {
            this.genomeID = genome_id;
            this.name = genome_name;
            this.distance = dist;
        }

        @Override
        public int compareTo(Neighbor o) {
            int retVal = Double.compare(this.distance, o.distance);
            if (retVal == 0)
                retVal = this.genomeID.compareTo(o.genomeID);
            return retVal;
        }

        /**
         * @return the ID of the neighbor genome
         */
        public String getID() {
            return this.genomeID;
        }

        /**
         * @return the name of the neighbor genome
         */
        public String getName() {
            return this.name;
        }

        /**
         * @return the distance to the representative
         */
        public double getDistance() {
            return this.distance;
        }

    }

    /**
     * This is a comparator used to sort the representative genomes from the largest neighborhoods to the smallest.
     */
    protected class RepSorter implements Comparator<String> {

        @Override
        public int compare(String o1, String o2) {
            // Get the neighbor lists.
            List<Neighbor> list1 = SyntheticSampleProcessor.this.neighborhoodMap.getOrDefault(o1, NO_NEIGHBORS);
            List<Neighbor> list2 = SyntheticSampleProcessor.this.neighborhoodMap.getOrDefault(o2, NO_NEIGHBORS);
            // Compare the lengths.  Longer sorts first.
            int retVal = list2.size() - list1.size();
            // Fall back to lexical by genome ID.
            if (retVal == 0)
                retVal = o1.compareTo(o2);
            return retVal;
        }

    }

    @Override
    protected void setReaderDefaults() {
        this.outFile = null;
        this.minGenomes = 2000;
        this.contigFrac = 1.0;
    }

    @Override
    protected void validateReaderParms() throws IOException, ParseFailureException {
        if (this.minGenomes < 0)
            throw new ParseFailureException("Minimum genome count cannot be negative.");
        if (this.contigFrac <= 0 || this.contigFrac > 1.0)
            throw new ParseFailureException("Contig fraction must be greater than 0 and less than or equal to 1.");
        // Connect to PATRIC.
        this.p3 = new P3Connection();
    }

    @Override
    protected void validateReaderInput(TabbedLineReader reader) throws IOException {
        // Verify we have all the fields we need.
        this.genomeIdColIdx = reader.findField("genome_id");
        this.genomeNameColIdx = reader.findField("genome_name");
        this.repIdColIdx = reader.findField("rep_id");
        this.distanceColIdx = reader.findField("distance");
    }

    @Override
    protected void runReader(TabbedLineReader reader) throws Exception {
        // Start by opening the output file.  We do this first so that if we can't write to it, we find out
        // before going through all the work of finding the contigs.
        try (FastaOutputStream outStream = this.getOutputStream()) {
            // Get some counters.
            int genomesIn = 0;
            int neighborsFound = 0;
            // Initialize the neighbor map.
            this.neighborhoodMap = new HashMap<String, List<Neighbor>>(2000);
            // Now we run through the list file, creating the neighbor lists.
            log.info("Reading input file.");
            for (var line : reader) {
                genomesIn++;
                // Get this genome's information.
                String genomeID = line.get(this.genomeIdColIdx);
                String repID = line.get(this.repIdColIdx);
                // Only proceed if the genome is not representing itself.
                if (! genomeID.contentEquals(repID)) {
                    // Get the remaining data and create a neighbor object.
                    String name = line.get(this.genomeNameColIdx);
                    double distance = line.getDouble(this.distanceColIdx);
                    var neighborInfo = new Neighbor(genomeID, name, distance);
                    var neighborList = this.neighborhoodMap.computeIfAbsent(repID, x -> new ArrayList<Neighbor>());
                    neighborList.add(neighborInfo);
                    neighborsFound++;
                }
            }
            log.info("{} genomes read, {} stored in neighbor map, {} representatives have neighbors.", genomesIn, neighborsFound,
                    this.neighborhoodMap.size());
            // Now we have the neighborhood map.  We must sort all the neighbor lists.  We will try to get genomes out of each
            // list at equally-spaced intervals to get a mix of distances.
            this.neighborhoodMap.values().forEach(x -> Collections.sort(x));
            log.info("Done sorting neighborhoods.");
            // Now we create a list of representative genome IDs.  We run through the list, counting the number of neighbors
            // we need to fulfill the minimum.  On each pass, we add one to the neighbor count of all remaining representatives.
            // If the count reaches the size of the neighborhood, we remove the representative from the list.  At the end
            // of a pass, we are done if the list is empty or we have reached the minimum.
            CountMap<String> proposedCounts = new CountMap<String>();
            int totalProposed = 0;
            var repList = new ArrayList<String>(this.neighborhoodMap.keySet());
            log.info("Scanning neighborhoods to plan genome selection.");
            while (repList.size() > 0 && totalProposed < this.minGenomes) {
                var iter = repList.iterator();
                while (iter.hasNext()) {
                    var genomeID = iter.next();
                    int counted = proposedCounts.count(genomeID);
                    totalProposed++;
                    // Remove this genome if we've used its whole neighborhood.
                    if (counted >= this.neighborhoodMap.get(genomeID).size())
                        iter.remove();
                }
            }
            log.info("{} genomes will be output.", totalProposed);
            // Now we run through the neighborhoods, selecting genomes.
            int gCount = 0;
            int seqCount = 0;
            for (var counter : proposedCounts.counts()) {
                String repID = counter.getKey();
                int count = counter.getCount();
                var neighborhood = this.neighborhoodMap.get(repID);
                Iterator<Neighbor> iter = new Spacer<Neighbor>(neighborhood, count);
                while (iter.hasNext()) {
                    var neighbor = iter.next();
                    var genomeID = neighbor.getID();
                    log.info("Processing genome {} for rep {}: {}", genomeID, repID, neighbor.getName());
                    gCount++;
                    // Get the contigs from PATRIC.
                    var gto = P3Genome.load(this.p3, repID, P3Genome.Details.CONTIGS);
                    for (Contig contig : gto.getContigs()) {
                        // Verify that we want this contig.
                        if (Math.random() < this.contigFrac) {
                            // Format the label and comment.
                            String label = genomeID + ":" + contig.getId();
                            String comment = neighbor.getName() + "\t" + repID + "\t" + Double.toString(neighbor.getDistance());
                            Sequence seq = new Sequence(label, comment, contig.getSequence());
                            outStream.write(seq);
                            seqCount++;
                        }
                    }
                }
            }
            log.info("{} genomes and {} sequences output.", gCount, seqCount);
        }
    }

    /**
     * @return the FASTA output stream to use for writing the contigs
     *
     * @throws FileNotFoundException
     */
    private FastaOutputStream getOutputStream() throws FileNotFoundException {
        FastaOutputStream retVal;
        if (this.outFile != null)
            retVal = new FastaOutputStream(this.outFile);
        else
            retVal = new FastaOutputStream(System.out);
        return retVal;
    }

}
