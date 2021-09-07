package com.hartwig.hmftools.linx.fusion;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.linx.fusion.rna.RnaFusionData;
import com.hartwig.hmftools.linx.germline.GermlinePonCache;

import org.apache.commons.cli.CommandLine;

public class FusionResources
{
    private final KnownFusionCache mKnownFusionCache;
    private final GermlinePonCache mGermlinePonCache;

    private final Map<String, List<RnaFusionData>> mSampleRnaData;

    public FusionResources(final CommandLine cmd)
    {
        mKnownFusionCache = new KnownFusionCache();
        mKnownFusionCache.loadFromFile(cmd);

        mGermlinePonCache = new GermlinePonCache(cmd);

        mSampleRnaData = Maps.newHashMap();
    }

    public KnownFusionCache knownFusionCache() { return mKnownFusionCache; }

    public GermlinePonCache germlinePonCache() { return mGermlinePonCache; }

    public final Map<String, List<RnaFusionData>> sampleRnaData() { return mSampleRnaData; }
}
