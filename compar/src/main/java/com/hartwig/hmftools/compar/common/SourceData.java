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

    public SourceData(final SourceType type, final PipelineSourcePaths pipelinePaths)
    {
        Type = type;
        PipelinePaths = pipelinePaths;
        Truthset = null;

        SampleIdMapping = Maps.newHashMap();
        ReferenceSampleIdMapping = Maps.newHashMap();
    }

    public SourceData(final SourceType type, final TruthsetCache truthset)
    {
        Type = type;
        PipelinePaths = null;
        Truthset = truthset;

        SampleIdMapping = Maps.newHashMap();
        ReferenceSampleIdMapping = Maps.newHashMap();
    }

    public String configName() { return Type.toString().toLowerCase(); }

    public String toString() { return Type.toString(); }
}
