package com.hartwig.hmftools.cup.cli;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.cup.prep.CuppaDataPrep;

public class PrepPlusPredictionMain {
    private static String[] extractArgs(String prefix, String[] allArgs) {
        List<String> allAsString = Arrays.stream(allArgs).collect(Collectors.toList());
        List<String> matchingArgs = new ArrayList<>();
        List<String> accumulator = new ArrayList<>();
        String argPrefix = String.format("-%s_", prefix);
        for (String arg: allAsString) {
            if (arg.startsWith(argPrefix)) {
                matchingArgs.addAll(accumulator);
                accumulator.clear();
                accumulator.add("-" + arg.substring(argPrefix.length()));
            }
            else {
                if (arg.startsWith("-")) {
                    matchingArgs.addAll(accumulator);
                    accumulator.clear();
                } else {
                    if (!accumulator.isEmpty()) {
                        accumulator.add(arg);
                    }
                }
            }
        }
        matchingArgs.addAll(accumulator);
        return matchingArgs.toArray(new String[] {});
    }

    public static void main(String[] args) {
        try {
            CuppaDataPrep.main(extractArgs("prep", args));
            PredictionRunner.main(extractArgs("prediction", args));
            System.exit(0);
        } catch (Exception e) {
            System.err.println("Failed to run combined Cuppa applications");
            e.printStackTrace();
            System.exit(1);
        }
    }
}
