package com.hartwig.hmftools.finding.datamodel;

import java.util.SortedMap;
import java.util.SortedSet;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

// there are two modes of prediction:
// CHORD is non cancer type specific (WGS),
// vCHORD is cancer type specific (panels)
@RecordBuilder
public record HomologousRecombination(
        @NotNull SortedMap<HrdCancerType, Prediction> predictions,
        @NotNull SortedSet<String> drivingGenes
)
{
    public enum HrdCancerType
    {
        PAN_CANCER,
        BREAST,
        OVARIAN,
        PANCREATIC,
        PROSTATE,
        OTHER
    }

    public enum Status
    {
        UNDETERMINED,
        HR_PROFICIENT,
        HR_DEFICIENT
    }

    public enum HrdType
    {
        BRACA1_TYPE,
        BRACA2_TYPE,
        CANNOT_BE_DETERMINED
    }

    @RecordBuilder
    public record Prediction(
            @NotNull String findingKey,
            @NotNull HrdCancerType cancerType,
            @NotNull Status status,
            @NotNull ThresholdValue hrdProbability,
            @Nullable Double brca1Probability,
            @Nullable Double brca2Probability,
            @Nullable HrdType hrdType)
            implements Finding
    {
    }
}
