package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.jetbrains.annotations.NotNull;

final class ReportingCopyNumberFilters {

    private static final double REL_GAIN = 3;
    private static final double ABS_LOSS = 0.5;

    private ReportingCopyNumberFilters() {
    }

    @NotNull
    static List<GeneCopyNumber> filterForReporting(@NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull GeneModel panelGeneModel) {
        return geneCopyNumbers.stream()
                .filter(copyNumber -> includeInReport(copyNumber.value(),
                        panelGeneModel.isAmplificationReportable(copyNumber.gene()),
                        panelGeneModel.isDeletionReportable(copyNumber.gene())))
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    static boolean includeInReport(double copyNumber, boolean isAmplificationReportable, boolean isDeletionReportable) {
        // KODU: Assume we only have significant events here.
        return (Doubles.greaterThan(copyNumber, ABS_LOSS) && isAmplificationReportable) ||
                (Doubles.lessOrEqual(copyNumber, ABS_LOSS) && isDeletionReportable);

    }
}
