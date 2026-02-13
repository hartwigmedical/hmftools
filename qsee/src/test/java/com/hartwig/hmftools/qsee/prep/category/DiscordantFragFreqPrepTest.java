package com.hartwig.hmftools.qsee.prep.category;


import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.common.sv.DiscordantFragType;
import com.hartwig.hmftools.common.sv.EsveeDiscordantStats;
import com.hartwig.hmftools.qsee.prep.category.discordant.DiscordantFragGroup;

import org.junit.Test;

public class DiscordantFragFreqPrepTest
{
    @Test
    public void canCalcDiscordantProportions()
    {
        long[] typeCounts = new long[DiscordantFragType.values().length];

        typeCounts[DiscordantFragType.Del1To5K.ordinal()] = 5;
        typeCounts[DiscordantFragType.Dup1To5K.ordinal()] = 5;

        EsveeDiscordantStats discordantStats = new EsveeDiscordantStats(100, Integer.MIN_VALUE, typeCounts);
        Map<DiscordantFragGroup, Double> discPropPerGroup = DiscordantFragFreqPrep.calcDiscordantProportions(discordantStats);

        assertEquals(0.1, discPropPerGroup.get(DiscordantFragGroup.DEL_DUP_SHORT), 0.001);
    }
}
