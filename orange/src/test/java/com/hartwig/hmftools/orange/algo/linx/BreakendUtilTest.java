package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BreakendUtilTest
{
    @Test
    public void canCreatePairsPerSv()
    {
        LinxBreakend breakend1 = LinxOrangeTestFactory.breakendBuilder().svId(1).build();
        LinxBreakend breakend2 = LinxOrangeTestFactory.breakendBuilder().svId(1).build();
        LinxBreakend breakend3 = LinxOrangeTestFactory.breakendBuilder().svId(2).build();
        LinxBreakend breakend4 = LinxOrangeTestFactory.breakendBuilder().svId(2).build();
        LinxBreakend breakend5 = LinxOrangeTestFactory.breakendBuilder().svId(3).build();

        List<LinxBreakend> breakends = Lists.newArrayList(breakend1, breakend2, breakend3, breakend4, breakend5);
        List<Pair<LinxBreakend, LinxBreakend>> pairs = BreakendUtil.createPairsPerSvId(breakends);

        assertEquals(2, pairs.size());
        Pair<LinxBreakend, LinxBreakend> pair1 = findBySvId(pairs, 1);
        assertTrue(pair1.getLeft().equals(breakend1) || pair1.getLeft().equals(breakend2));
        assertTrue(pair1.getRight().equals(breakend1) || pair1.getRight().equals(breakend2));

        Pair<LinxBreakend, LinxBreakend> pair2 = findBySvId(pairs, 2);
        assertTrue(pair2.getLeft().equals(breakend3) || pair2.getLeft().equals(breakend4));
        assertTrue(pair2.getRight().equals(breakend3) || pair2.getRight().equals(breakend4));
    }

    @NotNull
    private static Pair<LinxBreakend, LinxBreakend> findBySvId(@NotNull List<Pair<LinxBreakend, LinxBreakend>> pairs, int svIdToFind)
    {
        for(Pair<LinxBreakend, LinxBreakend> pair : pairs)
        {
            if(pair.getLeft().svId() == svIdToFind && pair.getRight().svId() == svIdToFind)
            {
                return pair;
            }
        }

        throw new IllegalStateException("Could not find breakend pair with svId: " + svIdToFind);
    }
}