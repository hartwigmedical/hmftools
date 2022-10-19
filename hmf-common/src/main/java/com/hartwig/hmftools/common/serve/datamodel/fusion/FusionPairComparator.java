package com.hartwig.hmftools.common.serve.datamodel.fusion;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FusionPairComparator implements Comparator<FusionPair> {

    @Override
    public int compare(@NotNull FusionPair pair1, @NotNull FusionPair pair2) {
        int geneUpCompare = pair1.geneUp().compareTo(pair2.geneUp());
        if (geneUpCompare != 0) {
            return geneUpCompare;
        }

        int geneDownCompare = pair1.geneDown().compareTo(pair2.geneDown());
        if (geneDownCompare != 0) {
            return geneDownCompare;
        }

        int exonUpCompare = compare(pair1.minExonUp(), pair2.minExonUp());
        if (exonUpCompare != 0) {
            return exonUpCompare;
        }

        return compare(pair1.minExonDown(), pair2.minExonDown());
    }

    private static int compare(@Nullable Integer int1, @Nullable Integer int2) {
        if (int1 == null && int2 == null) {
            return 0;
        } else if (int1 == null) {
            return 1;
        } else if (int2 == null) {
            return -1;
        }

        return int1.compareTo(int2);
    }
}
