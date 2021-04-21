package com.hartwig.hmftools.serve.extraction.exon;

import java.util.Comparator;

import com.hartwig.hmftools.serve.extraction.range.RangeCompare;

import org.jetbrains.annotations.NotNull;

class KnownExonComparator implements Comparator<KnownExon> {

    @Override
    public int compare(@NotNull KnownExon exon1, @NotNull KnownExon exon2) {
        int rangeCompare = RangeCompare.compare(exon1.annotation(), exon2.annotation());

        if (rangeCompare != 0) {
            return rangeCompare;
        }

        int transcriptCompare = exon1.annotation().transcript().compareTo(exon2.annotation().transcript());
        if (transcriptCompare != 0) {
            return transcriptCompare;
        }

        return Integer.compare(exon1.annotation().exonIndex(), exon2.annotation().exonIndex());
    }
}
