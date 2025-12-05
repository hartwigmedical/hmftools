package com.hartwig.hmftools.panelbuilder;

import org.jetbrains.annotations.Nullable;

// Metadata information about a region that we want to cover with probes.
public record TargetMetadata(
        Type type,
        String extraInfo,
        // TODO? make this a proper interface
        @Nullable Object extraData
)
{
    public TargetMetadata(final Type type, final String extraInfo)
    {
        this(type, extraInfo, null);
    }

    public enum Type
    {
        GENE,
        CN_BACKBONE,
        CDR3,
        CUSTOM_REGION,
        CUSTOM_SV,
        SAMPLE_SNV_INDEL_DRIVER,
        SAMPLE_SNV_INDEL_OTHER,
        SAMPLE_SV_AMP_DRIVER,
        SAMPLE_SV_DEL_DRIVER,
        SAMPLE_SV_FUSION_DRIVER,
        SAMPLE_SV_DISRUPTION_DRIVER,
        SAMPLE_GERMLINE_SNV_INDEL_DRIVER,
        SAMPLE_GERMLINE_SV_DRIVER,
    }
}
