package com.hartwig.hmftools.finding.datamodel;

import java.util.Set;

import jakarta.validation.constraints.NotNull;

@SuppressWarnings("unused")
@RecordBuilder
public record PurityPloidyFit(
        @NotNull Qc qc,
        @NotNull FittedPurityMethod fittedPurityMethod,
        double purity,
        double minPurity,
        double maxPurity,
        double ploidy,
        double minPloidy,
        double maxPloidy
)
{
    @RecordBuilder
    public record Qc(
            @NotNull Set<QCStatus> status,
            @NotNull Set<GermlineAberration> germlineAberrations,
            int amberMeanDepth,
            double contamination,
            int totalCopyNumberSegments,
            int unsupportedCopyNumberSegments,
            int deletedGenes
    )
    {
    }

    public enum FittedPurityMethod
    {
        NORMAL,
        HIGHLY_DIPLOID,
        SOMATIC,
        NO_TUMOR
    }

    public enum QCStatus
    {
        PASS,

        WARN_DELETED_GENES,
        WARN_HIGH_COPY_NUMBER_NOISE,
        WARN_GENDER_MISMATCH,
        WARN_LOW_PURITY,
        WARN_TINC,

        FAIL_CONTAMINATION,
        FAIL_NO_TUMOR,
        FAIL_TINC
    }

    public enum GermlineAberration
    {
        NONE,
        MOSAIC_X,
        KLINEFELTER,
        XYY,
        TRISOMY_X,
        TRISOMY_13,
        TRISOMY_15,
        TRISOMY_18,
        TRISOMY_21
    }

    public boolean isFail()
    {
        return isFailNoTumor() || isContaminated();
    }

    public boolean isContaminated()
    {
        return qc().status().contains(QCStatus.FAIL_CONTAMINATION);
    }

    public boolean isFailNoTumor()
    {
        return qc().status().contains(QCStatus.FAIL_NO_TUMOR);
    }

    public boolean containsTumorCells()
    {
        return !isFailNoTumor();
    }
}
