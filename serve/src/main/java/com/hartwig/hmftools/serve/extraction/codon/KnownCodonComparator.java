package com.hartwig.hmftools.serve.extraction.codon;

import java.util.Comparator;

import com.hartwig.hmftools.serve.extraction.range.RangeCompare;

import org.jetbrains.annotations.NotNull;

class KnownCodonComparator implements Comparator<KnownCodon> {

    @Override
    public int compare(@NotNull KnownCodon codon1, @NotNull KnownCodon codon2) {
        int rangeCompare = RangeCompare.compare(codon1.annotation(), codon2.annotation());

        if (rangeCompare != 0) {
            return rangeCompare;
        }

        int transcriptCompare = codon1.annotation().transcript().compareTo(codon2.annotation().transcript());
        if (transcriptCompare != 0) {
            return transcriptCompare;
        }

        return Integer.compare(codon1.annotation().codonIndex(), codon2.annotation().codonIndex());
    }
}
