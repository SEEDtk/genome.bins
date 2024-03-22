/**
 *
 */
package org.theseed.genome;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.theseed.stats.Shuffler;

/**
 * This object performs analysis of a source for synthetic sample genomes.  When a new source is
 * created, we just build a new subclass of this method.
 *
 * @author Bruce Parrello
 *
 */
public abstract class SynthFactory {

    // FIELDS
    /** maximum number of genomes to return */
    private int maxGenomes;
    /** number of genomes returned */
    private int nGenomes;


    /**
     * Construct the synthetic analyzer.
     *
     * @param max		maximum number of genomes to return
     */
    public SynthFactory(int max) {
        this.maxGenomes = max;
        this.nGenomes = 0;
    }

    /**
     * @return the next eligible genome, or NULL if there are none left
     *
     * @throws IOException
     */
    public Genome getNext() throws IOException {
        Genome retVal;
        if (this.nGenomes >= this.maxGenomes)
            retVal = null;
        else {
            retVal = this.getNextGenome();
            this.nGenomes++;
        }
        return retVal;
    }

    /**
     * @return the next eligible genome
     *
     * @throws IOException
     */
    protected abstract Genome getNextGenome() throws IOException;

    /**
     * Select a random set of genomes from all possibilities to insure we are under the maximum.
     * This is a utility method for synthetic sources that can do fast pre-selection.
     *
     * @param genomes		full set of genome descriptors
     *
     * @return a smaller list of genome descriptors to use
     */
    protected <T> List<T> selectGenomes(Collection<T> genomes) {
        // Randomly shuffle the genome ID collection to select the desired number of genomes.
        var buffer = new Shuffler<T>(genomes);
        buffer.shuffle(this.maxGenomes);
        // If the list is already small enough, we are done; otherwise, trim off the start of the list.
        List<T> retVal;
        if (buffer.size() <= this.maxGenomes)
            retVal = buffer;
        else
            retVal = new ArrayList<T>(buffer.subList(0, this.maxGenomes));
        return retVal;
    }

}
