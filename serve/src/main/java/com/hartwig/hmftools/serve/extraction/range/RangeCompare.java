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

        int transcriptCompare = annotation1.transcript().compareTo(annotation2.transcript());
        if (transcriptCompare != 0) {
            return transcriptCompare;
        }

        int rankCompare = Integer.compare(annotation1.rank(), annotation2.rank());
        if (rankCompare != 0) {
            return rankCompare;
        }

        return annotation1.mutationType().compareTo(annotation2.mutationType());
    }
}
