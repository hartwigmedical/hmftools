package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.actionability.cnv.SignificantGeneCopyNumberFilter;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.jetbrains.annotations.NotNull;

final class ReportingCopyNumberFilters {

    private ReportingCopyNumberFilters() {
    }

    @NotNull
    static List<GeneCopyNumber> filterForReporting(@NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull GeneModel panelGeneModel,
            @NotNull Gender gender, double averageTumorPloidy) {
        List<GeneCopyNumber> significantGeneCopyNumbers =
                SignificantGeneCopyNumberFilter.filterForSignificance(geneCopyNumbers, averageTumorPloidy);

        return significantGeneCopyNumbers.stream()
                .filter(copyNumber -> includeInReport(copyNumber.minCopyNumber(),
                        HumanChromosome.valueOf(copyNumber).isDiploid(gender),
                        panelGeneModel.isAmplificationReportable(copyNumber.gene()),
                        panelGeneModel.isDeletionReportable(copyNumber.gene())))
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    static boolean includeInReport(double significantCopyNumber, boolean isDiploidChromosome, boolean isAmplificationReportable,
            boolean isDeletionReportable) {
        // Assume we only have significant events here.
        double normalPloidy = isDiploidChromosome ? 2D : 1D;
        return (significantCopyNumber > normalPloidy && isAmplificationReportable) || (significantCopyNumber < normalPloidy
                && isDeletionReportable);
    }
}
