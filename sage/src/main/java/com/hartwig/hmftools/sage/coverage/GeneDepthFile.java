package com.hartwig.hmftools.sage.coverage;

import java.util.StringJoiner;

import org.jetbrains.annotations.NotNull;

public class GeneDepthFile {
    private static final String DELIMITER = "\t";

    static String header(@NotNull final GeneDepth depth) {
        StringJoiner joiner = new StringJoiner(DELIMITER).add("gene");
        int[] depthCounts = depth.depthCounts();
        for (int i = 0; i < depthCounts.length; i++) {
            if (i == depthCounts.length - 1) {
                joiner.add(i + "+");
            } else {
                joiner.add(String.valueOf(i));
            }
        }

        return joiner.toString();
    }

    @NotNull
    static String toString(@NotNull final GeneDepth depth) {
        StringJoiner joiner = new StringJoiner(DELIMITER).add(depth.gene());
        for (int i : depth.depthCounts()) {
            joiner.add(String.valueOf(i));
        }

        return joiner.toString();
    }

}
