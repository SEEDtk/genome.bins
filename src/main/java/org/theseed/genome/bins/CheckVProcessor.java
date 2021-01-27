/**
 *
 */
package org.theseed.genome.bins;

import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command processes the CheckV output for a binning run.  It will sort contigs into virus bins (vbin.fa) and produce a tab-delimited
 * report on the results.
 *
 * The positional parameters are the name of the directory containing the CheckV database, the name of the binning output directory,
 * and the name of the directory to contain the output files from this command.
 *
 * In the binning output directory
 *
 * 		output.contigs2reads.txt	contains the coverage for each contig ID in a 2-column, tab-delimited form.
 * 		unbinned.fasta				contains the actual contigs (FASTA)
 * 		CheckV/completeness.tsv		contains the CheckV output (tab-delimited)
 *
 * In the CheckV database, the genome_db directory contains a file "checkv_reps.tsv" that maps each
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -e	maximum error value for a contig to be considered for inclusion in a bin
 *
 *
 * @author Bruce Parrello
 *
 */
public class CheckVProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(CheckVProcessor.class);


    @Override
    protected void setDefaults() {
        //  code for setDefaults

    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        //  code for validateParms
        return false;
    }

    @Override
    protected void runCommand() throws Exception {
        //  code for runCommand

    }

}
