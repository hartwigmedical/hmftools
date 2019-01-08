package com.hartwig.hmftools.common.copynumber;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.jetbrains.annotations.NotNull;

public class FilterSignificantGeneCopyNumbers {

    private static final double REL_GAIN = 3;
    private static final double ABS_LOSS = 0.5;

    @NotNull
    public static List<GeneCopyNumber> filterForSignificance(@NotNull List<GeneCopyNumber> geneCopyNumbers, double averageTumorPloidy) {
        return geneCopyNumbers.stream()
                .filter(copyNumber -> isSignificant(averageTumorPloidy, copyNumber.value()))
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