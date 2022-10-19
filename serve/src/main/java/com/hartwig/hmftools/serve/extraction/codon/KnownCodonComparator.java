package com.hartwig.hmftools.serve.extraction.codon;

import java.util.Comparator;

import com.hartwig.hmftools.common.serve.datamodel.range.RangeCompare;

import org.jetbrains.annotations.NotNull;

class KnownCodonComparator implements Comparator<KnownCodon> {

    @Override
    public int compare(@NotNull KnownCodon codon1, @NotNull KnownCodon codon2) {
        return RangeCompare.compare(codon1.annotation(), codon2.annotation());
    }
}
