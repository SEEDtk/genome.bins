/**
 *
 */
package org.theseed.genome;

import java.io.File;
import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;

/**
 * This is utility object used in writing synthetic samples.  The constructor creates the FASTA output stream,
 * and the main output method writes the sample sequences for a particular genome.
 *
 * It is auto-closeable so that the file is closed properly.
 *
 * @author Bruce Parrello
 *
 */
public class SynthWriter implements AutoCloseable {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SynthWriter.class);
    /** output stream */
    private FastaOutputStream outStream;
    /** fraction of contigs to output */
    private double contigFrac;
    /** number of sequences written */
    private int outCount;
    /** number of genomes written */
    private int gCount;

    /**
     * Construct a synthetic sample writer.
     *
     * @param outFile	output file, or NULL to use the standard output
     * @param frac		fraction (from 0 to 1) of contigs to output
     *
     * @throws IOException
     */
    public SynthWriter(File outFile, double frac) throws IOException {
        this.contigFrac = frac;
        if (outFile == null) {
            log.info("Sequences will be written to the standard output.");
            this.outStream = new FastaOutputStream(System.out);
        } else {
            log.info("Sequences will be written to {}.", outFile);
            this.outStream = new FastaOutputStream(outFile);
        }
        this.outCount = 0;
        this.gCount = 0;
    }

    /**
     * Write the contigs from a genome that is a specific distance from a representative.
     *
     * @param	gto		gto to write
     * @param	repId	representative genome's ID
     * @param	dist	distance
     *
     * @throws IOException
     */
    public void writeGenome(Genome gto, String repId, double dist) throws IOException {
        String genomeID = gto.getId();
        String name = gto.getName();
        log.info("Writing contigs from genome {}.", gto);
        for (Contig contig : gto.getContigs()) {
            // Verify that we want this contig.
            if (Math.random() < this.contigFrac) {
                // Format the label and comment.
                String label = genomeID + ":" + contig.getId();
                String comment = name + "\t" + repId + "\t" + Double.toString(dist);
                Sequence seq = new Sequence(label, comment, contig.getSequence());
                this.outStream.write(seq);
                this.outCount++;
            }
        }
        this.gCount++;
    }

    @Override
    public void close() {
        this.outStream.close();
        log.info("{} total sequences written from {} genomes.", this.outCount, this.gCount);
    }


}
