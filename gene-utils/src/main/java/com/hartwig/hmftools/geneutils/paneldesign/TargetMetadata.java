package com.hartwig.hmftools.geneutils.paneldesign;

import org.jetbrains.annotations.Nullable;

// Metadata information about a region that we want to cover with probes.
public record TargetMetadata(
        Type type,
        String extraInfo,
        @Nullable ExtraData extraData
)
{
    public TargetMetadata(Type type, String extraInfo)
    {
        this(type, extraInfo, null);
    }

    public enum Type
    {
        GENE,
        CN_BACKBONE,
        CUSTOM
    }

    public sealed interface ExtraData permits CopyNumberBackbone.TargetMetadataExtra
    {
    }
}
