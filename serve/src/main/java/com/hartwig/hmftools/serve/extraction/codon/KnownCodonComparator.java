package com.hartwig.hmftools.serve.extraction.codon;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

class KnownCodonComparator implements Comparator<KnownCodon> {

    @Override
    public int compare(@NotNull KnownCodon codon1, @NotNull KnownCodon codon2) {
        return codon1.annotation().compareTo(codon2.annotation());
    }
}
