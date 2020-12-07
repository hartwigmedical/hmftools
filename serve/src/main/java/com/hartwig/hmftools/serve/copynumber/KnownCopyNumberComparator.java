package com.hartwig.hmftools.serve.copynumber;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

class KnownCopyNumberComparator implements Comparator<KnownCopyNumber> {

    @Override
    public int compare(@NotNull KnownCopyNumber copyNumber1, @NotNull KnownCopyNumber copyNumber2) {
        int geneCompare = copyNumber1.gene().compareTo(copyNumber2.gene());
        if (geneCompare != 0) {
            return geneCompare;
        }

        return copyNumber1.type().toString().compareTo(copyNumber2.type().toString());
    }
}
