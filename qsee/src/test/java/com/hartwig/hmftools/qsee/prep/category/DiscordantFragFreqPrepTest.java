package com.hartwig.hmftools.qsee.prep.category;

import static org.junit.Assert.assertEquals;

import java.util.EnumMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.sv.DiscordantFragType;
import com.hartwig.hmftools.common.sv.EsveeDiscordantStats;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.category.discordant.DiscordantFragGroup;

import org.junit.Test;

public class DiscordantFragFreqPrepTest
{
    private static final int DUMMY_READ_COUNT = 1;
    private static final int DUMMY_TOTAL_READS = 100;
    private static final double DUMMY_BASE_PROPORTION = DUMMY_READ_COUNT / (double) DUMMY_TOTAL_READS;

    private static final Map<DiscordantFragGroup, Integer> TYPES_PER_GROUP = getTypesPerGroup();

    @Test
    public void canCalcDiscordantProportions()
    {
        long[] typeCounts = new long[DiscordantFragType.values().length];
        for(DiscordantFragType type : DiscordantFragType.values())
        {
            typeCounts[type.ordinal()] = DUMMY_READ_COUNT;
        }

        EsveeDiscordantStats discordantStats = new EsveeDiscordantStats(DUMMY_TOTAL_READS, Integer.MIN_VALUE, typeCounts);

        Map<DiscordantFragGroup, Double> discPropPerGroup = DiscordantFragFreqPrep.calcDiscordantProportions(discordantStats);
        List<Feature> actualFeatures = DiscordantFragFreqPrep.formFeatures(discPropPerGroup);

        List<Feature> expectedFeatures = List.of(
                createExpectedFeature(DiscordantFragGroup.DEL_DUP_SHORT),
                createExpectedFeature(DiscordantFragGroup.DEL_DUP_MEDIUM),
                createExpectedFeature(DiscordantFragGroup.DEL_DUP_LONG),
                createExpectedFeature(DiscordantFragGroup.INV_SHORT),
                createExpectedFeature(DiscordantFragGroup.INV_MEDIUM),
                createExpectedFeature(DiscordantFragGroup.INV_LONG),
                createExpectedFeature(DiscordantFragGroup.TRANSLOCATION)
        );

        for(int i = 0; i < expectedFeatures.size(); i++)
        {
            Feature actualFeature = actualFeatures.get(i);
            Feature expectedFeature = expectedFeatures.get(i);

            assertEquals(expectedFeature.key().name(), actualFeature.key().name());
            assertEquals(expectedFeature.value(), actualFeature.value(), 0.001);
        }
    }

    private static Feature createExpectedFeature(DiscordantFragGroup discordantFragGroup)
    {
        String featureName = "DiscordantFragType=" + discordantFragGroup.getName();
        double featureValue = TYPES_PER_GROUP.get(discordantFragGroup) * DUMMY_BASE_PROPORTION;

        return new Feature(featureName, featureValue, FeatureType.DISCORDANT_FRAG_FREQ, SourceTool.ESVEE, null);
    }

    private static Map<DiscordantFragGroup, Integer> getTypesPerGroup()
    {
        Map<DiscordantFragGroup, Integer> typeCountsPerGroup = new EnumMap<>(DiscordantFragGroup.class);

        for(DiscordantFragType type : DiscordantFragType.values())
        {
            DiscordantFragGroup group = DiscordantFragGroup.fromType(type);

            if(!typeCountsPerGroup.containsKey(group))
                typeCountsPerGroup.put(group, 0);

            typeCountsPerGroup.put(group, typeCountsPerGroup.get(group) + 1);
        }

        return typeCountsPerGroup;
    }
}
