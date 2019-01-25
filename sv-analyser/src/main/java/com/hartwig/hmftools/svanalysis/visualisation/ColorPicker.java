package com.hartwig.hmftools.svanalysis.visualisation;

import org.jetbrains.annotations.NotNull;

public interface ColorPicker {

    String[] COLOURS =
            new String[] { "166,206,227", "106,61,154", "51,160,44", "227,26,28", "255,127,0", "31,120,180", "202,178,214", "178,223,138",
                    "251,154,153", "253,191,111", "255,255,153" };

    @NotNull
    default String color(final int clusterId, final int chainId) {
        return chainId < COLOURS.length ? "color=" + COLOURS[chainId] : "color=black";
    }

}
