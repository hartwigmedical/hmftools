package com.hartwig.hmftools.qsee.prep.category;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.sv.DiscordantFragType;
import com.hartwig.hmftools.common.sv.EsveeDiscordantStats;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;
import com.hartwig.hmftools.qsee.prep.discordant.DiscordantFragGroup;

import org.jetbrains.annotations.NotNull;

public class DiscordantFragFreqPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.ESVEE;

    private static final String FIELD_FRAG_TYPE = "DiscordantFragType";

    public DiscordantFragFreqPrep(CommonPrepConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }

    private EsveeDiscordantStats loadDiscordantStats(String sampleId)
    {
        String baseDir = mConfig.getEsveeDir(sampleId);
        String filePath = EsveeDiscordantStats.generateFilename(baseDir, sampleId);
        return EsveeDiscordantStats.read(filePath);
    }

    static Map<DiscordantFragGroup, Double> calcDiscordantProportions(EsveeDiscordantStats discordantStats)
    {
        Map<DiscordantFragGroup, Long> discCountPerGroup = new EnumMap<>(DiscordantFragGroup.class);
        for(int i = 0; i < discordantStats.TypeCounts.length; i++)
        {
            DiscordantFragType type = DiscordantFragType.values()[i];
            DiscordantFragGroup group = DiscordantFragGroup.fromType(type);
            long count = discordantStats.TypeCounts[i];

            long currentGroupCount = discCountPerGroup.getOrDefault(group, 0L);
            currentGroupCount += count;
            discCountPerGroup.put(group, currentGroupCount);
        }

        Map<DiscordantFragGroup, Double> discPropPerGroup = new EnumMap<>(DiscordantFragGroup.class);
        for(DiscordantFragGroup group : discCountPerGroup.keySet())
        {
            double discProp = (double) discCountPerGroup.get(group) / discordantStats.TotalReads;
            discPropPerGroup.put(group, discProp);
        }

        return discPropPerGroup;
    }

    private static List<Feature> formFeatures(Map<DiscordantFragGroup, Double> discPropPerGroup)
    {
        List<Feature> features = new ArrayList<>();
        for(DiscordantFragGroup group : discPropPerGroup.keySet())
        {
            String featureName = FeatureKey.formSingleFieldName(FIELD_FRAG_TYPE, group.getName());
            FeatureKey key = new FeatureKey(featureName, FeatureType.DISCORDANT_FRAG_FREQ, SOURCE_TOOL);
            Feature feature = new Feature(key, discPropPerGroup.get(group));
            features.add(feature);
        }

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        if(sampleType != SampleType.TUMOR)
        {
            return List.of();
        }

        EsveeDiscordantStats discordantStats = loadDiscordantStats(sampleId);
        Map<DiscordantFragGroup, Double> discPropPerGroup = calcDiscordantProportions(discordantStats);
        List<Feature> features = formFeatures(discPropPerGroup);

        return features;
    }

}
