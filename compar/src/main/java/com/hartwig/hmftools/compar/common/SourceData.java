package com.hartwig.hmftools.compar.common;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.Nullable;

public class SourceData
{
    public final SourceType Type;

    @Nullable
    public final PipelineSourcePaths PipelinePaths;

    @Nullable
    public final ExpectedCache Expected;

    @Nullable
    public final String Label;

    public final Map<String,String> SampleIdMapping; // maps from the common tumor SampleId to the tumor Sample Id for this source
    public final Map<String,String> ReferenceSampleIdMapping; // as above for the reference Id

    public SourceData(final SourceType type, final PipelineSourcePaths pipelinePaths, final ExpectedCache expected, final String label)
    {
        Type = type;
        PipelinePaths = pipelinePaths;
        Expected = expected;
        Label = label;

        SampleIdMapping = Maps.newHashMap();
        ReferenceSampleIdMapping = Maps.newHashMap();
    }

    public boolean isPipelineSourced() { return PipelinePaths != null; }
    public boolean isExpectedSourced() { return Expected != null; }

    public static SourceData fromPipelineSource(final SourceType type, PipelineSourcePaths pipelinePaths, final String label)
    {
        return new SourceData(type, pipelinePaths, null, label);
    }

    public static SourceData fromExpectedSource(final SourceType type, final ExpectedCache expected, final String label)
    {
        return new SourceData(type, null, expected, label);
    }

    public String configName() { return Type.toString().toLowerCase(); }

    public String toString() { return Type.toString(); }
}
