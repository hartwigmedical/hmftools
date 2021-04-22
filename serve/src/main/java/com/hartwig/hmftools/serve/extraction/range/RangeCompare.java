package com.hartwig.hmftools.serve.extraction.range;

import org.jetbrains.annotations.NotNull;

public final class RangeCompare {

    private RangeCompare() {
    }

    public static int compare(@NotNull RangeAnnotation annotation1, @NotNull RangeAnnotation annotation2) {
        int rangeCompare = annotation1.compareTo(annotation2);

        if (rangeCompare != 0) {
            return rangeCompare;
        }

        int geneCompare = annotation1.gene().compareTo(annotation2.gene());
        if (geneCompare != 0) {
            return geneCompare;
        }

        return annotation1.mutationType().compareTo(annotation2.mutationType());
    }
}
