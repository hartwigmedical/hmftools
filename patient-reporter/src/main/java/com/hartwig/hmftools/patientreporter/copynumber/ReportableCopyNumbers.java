package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

final class ReportableCopyNumbers {

    @VisibleForTesting
    static final double REL_GAIN = 3;
    @VisibleForTesting
    static final double ABS_LOSS = 0.5;

    private ReportableCopyNumbers() {
    }

    @NotNull
    static List<GeneCopyNumber> filterCopyNumbersForReport(final double samplePloidy, @NotNull final List<GeneCopyNumber> geneCopyNumbers) {
        return geneCopyNumbers.stream()
                .filter(copyNumber -> includeInReport(samplePloidy, copyNumber.value()))
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    static boolean includeInReport(final double samplePloidy, final double copyNumber) {
        if (Doubles.lessOrEqual(copyNumber, ABS_LOSS)) {
            return true;
        }

        double relativeCopyNumber = copyNumber / samplePloidy;
        return Doubles.greaterOrEqual(relativeCopyNumber, REL_GAIN);
    }
}
