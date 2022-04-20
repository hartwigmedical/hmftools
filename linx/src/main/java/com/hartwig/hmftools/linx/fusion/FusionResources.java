package com.hartwig.hmftools.linx.fusion;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.linx.fusion.rna.RnaFusionData;

import org.apache.commons.cli.CommandLine;

public class FusionResources
{
    private final KnownFusionCache mKnownFusionCache;

    private final Map<String, List<RnaFusionData>> mSampleRnaData;

    public FusionResources(final CommandLine cmd)
    {
        mKnownFusionCache = new KnownFusionCache();
        mKnownFusionCache.loadFromFile(cmd);

        mSampleRnaData = Maps.newHashMap();
    }

    public KnownFusionCache knownFusionCache() { return mKnownFusionCache; }

    public final Map<String, List<RnaFusionData>> sampleRnaData() { return mSampleRnaData; }
}
