package com.hartwig.hmftools.finding.datamodel;

import java.util.SortedSet;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@SuppressWarnings("unused")
@RecordBuilder
public record Qc(
        boolean isPass,
        @NotNull SortedSet<QCStatus> warnings,
        @NotNull SortedSet<QCStatus> errors,
        @NotNull SortedSet<GermlineAberration> germlineAberrations,
        int amberMeanDepth,
        double contamination,
        int totalCopyNumberSegments,
        int unsupportedCopyNumberSegments,
        int deletedGenes,
        @Nullable VisualisationFile visualisationFile
)
{
    public enum QCStatus
    {
        DELETED_GENES,
        HIGH_COPY_NUMBER_NOISE,
        GENDER_MISMATCH,
        LOW_PURITY,
        TUMOR_LOW_COVERAGE,
        TUMOR_LOW_MAPPED_PROPORTION,
        TUMOR_LOW_BASE_QUAL,
        TUMOR_LOW_MAP_QUAL,
        NORMAL_LOW_COVERAGE,
        NORMAL_LOW_MAPPED_PROPORTION,
        NORMAL_LOW_BASE_QUAL,
        NORMAL_LOW_MAP_QUAL,

        CONTAMINATION,
        NO_TUMOR,
        TUMOR_IN_NORMAL_CONTAMINATION
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
}
