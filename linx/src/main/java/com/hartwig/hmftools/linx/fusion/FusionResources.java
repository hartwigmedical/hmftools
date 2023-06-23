package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.linx.fusion.rna.RnaFusionData;

import org.apache.commons.cli.CommandLine;

public class FusionResources
{
    private final KnownFusionCache mKnownFusionCache;

    private final Map<String, List<RnaFusionData>> mSampleRnaData;

    public FusionResources(final ConfigBuilder configBuilder)
    {
        mKnownFusionCache = new KnownFusionCache();
        mKnownFusionCache.loadFromFile(configBuilder.getValue(KNOWN_FUSIONS_FILE));

        mSampleRnaData = Maps.newHashMap();
    }

    public KnownFusionCache knownFusionCache() { return mKnownFusionCache; }

    public final Map<String, List<RnaFusionData>> sampleRnaData() { return mSampleRnaData; }
}
