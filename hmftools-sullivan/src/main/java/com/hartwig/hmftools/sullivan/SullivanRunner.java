package com.hartwig.hmftools.sullivan;

public class SullivanRunner {

    public static void main(String[] args) {
        if (args.length != 2) {
            System.out.println("Usage <originalFastqPath> <recreatedFastqPath>");
        } else {
            SullivanAlgo.runSullivanAlgo(args[0], args[1]);
        }
    }
}
