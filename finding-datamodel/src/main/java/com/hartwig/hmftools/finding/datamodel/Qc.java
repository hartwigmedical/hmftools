package com.hartwig.hmftools.finding.datamodel;

import java.util.Set;

import jakarta.validation.constraints.NotNull;

@SuppressWarnings("unused")
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
    public enum QCStatus
    {
        PASS,

        WARN_DELETED_GENES,
        WARN_HIGH_COPY_NUMBER_NOISE,
        WARN_GENDER_MISMATCH,
        WARN_LOW_PURITY,
        WARN_TUMOR_IN_NORMAL_CONTAMINATION,
        WARN_TUMOR_LOW_COVERAGE,
        WARN_TUMOR_LOW_MAPPED_PROPORTION,
        WARN_TUMOR_LOW_BASE_QUAL,
        WARN_TUMOR_LOW_MAP_QUAL,
        WARN_NORMAL_LOW_COVERAGE,
        WARN_NORMAL_LOW_MAPPED_PROPORTION,
        WARN_NORMAL_LOW_BASE_QUAL,
        WARN_NORMAL_LOW_MAP_QUAL,

        FAIL_CONTAMINATION,
        FAIL_NO_TUMOR,
        FAIL_TUMOR_IN_NORMAL_CONTAMINATION
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
        return status().contains(QCStatus.FAIL_CONTAMINATION);
    }

    public boolean isLowPurity()
    {
        return status().contains(QCStatus.WARN_LOW_PURITY);
    }

    public boolean isFailNoTumor()
    {
        return status().contains(QCStatus.FAIL_NO_TUMOR);
    }

    public boolean containsTumorCells()
    {
        return !isFailNoTumor();
    }
}
