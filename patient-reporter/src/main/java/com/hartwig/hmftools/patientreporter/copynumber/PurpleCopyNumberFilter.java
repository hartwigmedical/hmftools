package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

final class PurpleCopyNumberFilter {

    @VisibleForTesting
    static final double ABS_GAIN = 8;
    @VisibleForTesting
    static final double REL_GAIN = 2.2;
    @VisibleForTesting
    static final double ABS_LOSS = 0.5;

    private PurpleCopyNumberFilter() {
    }

    @NotNull
    static List<GeneCopyNumber> filterCopyNumbersForReport(final double samplePloidy, @NotNull final List<GeneCopyNumber> geneCopyNumbers) {
        final List<GeneCopyNumber> result = Lists.newArrayList();
        for (final GeneCopyNumber geneCopyNumber : geneCopyNumbers) {
            if (includeInReport(samplePloidy, geneCopyNumber.minCopyNumber())) {
                result.add(ImmutableGeneCopyNumber.copyOf(geneCopyNumber));
            }
        }
        return result;
    }

    @VisibleForTesting
    static boolean includeInReport(final double samplePloidy, final double copyNumber) {
        if (Doubles.lessOrEqual(copyNumber, ABS_LOSS)) {
            return true;
        }

        if (Doubles.greaterOrEqual(copyNumber, ABS_GAIN)) {
            return true;
        }

        double relativeCopyNumber = copyNumber / samplePloidy;
        return Doubles.greaterOrEqual(relativeCopyNumber, REL_GAIN);
    }
}
