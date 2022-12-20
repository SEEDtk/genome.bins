package org.theseed.genome.bins;

import java.util.Arrays;

import org.theseed.utils.BaseProcessor;

/**
 * This performs various analyses of binning runs.  The possible commands are:
 *
 * bins		produce a report comparing bin quality to its distance to the reference genome
 * checkv	produce a report on the checkv output
 * repSynth	generate a synthetic sample from a repgen list file
 * synth	generate a synthetic sample from multiple sources
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        switch (command) {
        case "bins" :
            processor = new BinCheckProcessor();
            break;
        case "checkv" :
            processor = new CheckVProcessor();
            break;
        case "repSynth" :
            processor = new RepSyntheticSampleProcessor();
            break;
        case "synth" :
            processor = new SyntheticSampleProcessor();
            break;
        default :
            throw new RuntimeException("Invalid command " + command + ".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
