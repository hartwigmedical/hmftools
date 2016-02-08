package com.hartwig.hmftools.sullivan;

public class SullivanRunner {

    public static void main(String[] args) {
        if (args.length < 2 || args.length > 3) {
            System.out.println("Usage <originalFastqPath> <recreatedFastqPath> <optional:numRecordsToCheck>");
        } else {
            if (args.length == 3) {
                int numRecords = Integer.parseInt(args[2]);
                SullivanAlgo.runSullivanAlgo(args[0], args[1], numRecords);
            } else {
                SullivanAlgo.runSullivanAlgo(args[0], args[1]);
            }
        }
    }
}
