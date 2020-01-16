package com.hartwig.hmftools.protect.actionability.cnv;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public final class SignificantGeneCopyNumberFilter {

    private static final double REL_GAIN = 3;
    private static final double ABS_LOSS = 0.5;

    private SignificantGeneCopyNumberFilter() {
    }

    @NotNull
    public static List<GeneCopyNumber> filterForSignificance(@NotNull List<GeneCopyNumber> geneCopyNumbers, double averageTumorPloidy) {
        return geneCopyNumbers.stream()
                .filter(copyNumber -> isSignificant(averageTumorPloidy, copyNumber.minCopyNumber()))
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    static boolean isSignificant(double averageTumorPloidy, double copyNumber) {
        if (Doubles.lessOrEqual(copyNumber, ABS_LOSS)) {
            return true;
        }

        double relativeCopyNumber = copyNumber / averageTumorPloidy;
        return Doubles.greaterOrEqual(relativeCopyNumber, REL_GAIN);
    }
}