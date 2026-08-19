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
    public final TruthsetCache Truthset;

    public final Map<String,String> SampleIdMapping; // maps from the common tumor SampleId to the tumor Sample Id for this source
    public final Map<String,String> ReferenceSampleIdMapping; // as above for the reference Id

    public SourceData(final SourceType type, final PipelineSourcePaths pipelinePaths, final TruthsetCache truthset)
    {
        Type = type;
        PipelinePaths = pipelinePaths;
        Truthset = truthset;

        SampleIdMapping = Maps.newHashMap();
        ReferenceSampleIdMapping = Maps.newHashMap();
    }

    public static SourceData fromPipelineSource(final SourceType type, PipelineSourcePaths pipelinePaths)
    {
        return new SourceData(type, pipelinePaths, null);
    }

    public static SourceData fromTruthsetSource(final SourceType type, final TruthsetCache truthset)
    {
        return new SourceData(type, null, truthset);
    }

    public String configName() { return Type.toString().toLowerCase(); }

    public String toString() { return Type.toString(); }
}
