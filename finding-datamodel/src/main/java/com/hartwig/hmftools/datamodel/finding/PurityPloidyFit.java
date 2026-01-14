package com.hartwig.hmftools.datamodel.finding;

import java.util.Set;

import com.hartwig.hmftools.datamodel.purple.PurpleFittedPurityMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;

import org.jetbrains.annotations.NotNull;

public record PurityPloidyFit(
        @NotNull Qc qc,
        @NotNull PurpleFittedPurityMethod fittedPurityMethod,
        double purity,
        double minPurity,
        double maxPurity,
        double ploidy,
        double minPloidy,
        double maxPloidy
) {
    public record Qc(
        @NotNull Set<PurpleQCStatus> status,
        @NotNull Set<PurpleGermlineAberration> germlineAberrations,
        int amberMeanDepth,
        double contamination,
        int totalCopyNumberSegments,
        int unsupportedCopyNumberSegments,
        int deletedGenes
    ) {}

    public boolean containsTumorCells()
    {
        return !qc().status().contains(PurpleQCStatus.FAIL_NO_TUMOR);
    }
}
