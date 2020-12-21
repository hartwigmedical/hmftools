package com.hartwig.hmftools.serve.extraction.exon;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

class KnownExonComparator implements Comparator<KnownExon> {

    @Override
    public int compare(@NotNull KnownExon exon1, @NotNull KnownExon exon2) {
        return exon1.annotation().compareTo(exon2.annotation());
    }
}
