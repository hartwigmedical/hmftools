package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.calcEffectiveLength;

import static junit.framework.TestCase.assertEquals;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.common.BaseDepth;

import org.junit.Test;

public class MiscTest
{
    @Test
    public void testEffectiveLength()
    {
        final List<FragmentSize> fragmentLengthData = Lists.newArrayList();

        fragmentLengthData.add(new FragmentSize(100, 1));
        fragmentLengthData.add(new FragmentSize(200, 2));
        fragmentLengthData.add(new FragmentSize(300, 1));

        assertEquals(800, calcEffectiveLength(1000, fragmentLengthData), 0.001);

        assertEquals(83.33, calcEffectiveLength(250, fragmentLengthData), 0.1);
    }

    @Test
    public void testBaseDepth()
    {
        BaseDepth baseDepth = new BaseDepth();

        int[] baseRange = new int[SE_PAIR];

        baseRange[SE_START] = 100;
        baseRange[SE_END] = 200;

        baseDepth.initialise(baseRange);

        List<int[]> readCoords = Lists.newArrayList();
        readCoords.add(new int[] {100, 200});
        baseDepth.processRead(readCoords);

        readCoords.clear();
        readCoords.add(new int[] {100, 150});
        baseDepth.processRead(readCoords);

        // reads out the range have no effect
        readCoords.clear();
        readCoords.add(new int[] {50, 60});
        readCoords.add(new int[] {250, 260});
        baseDepth.processRead(readCoords);

        assertEquals(51, baseDepth.basesWithDepth());
        assertEquals(2, baseDepth.depthAtBase(100));
        assertEquals(2, baseDepth.depthAtBase(150));
        assertEquals(1, baseDepth.depthAtBase(200));
        assertEquals(0, baseDepth.depthAtBase(99));
        assertEquals(0, baseDepth.depthAtBase(201));

        readCoords.clear();
        readCoords.add(new int[] {100, 110});
        readCoords.add(new int[] {120, 130});
        baseDepth.processRead(readCoords);

        assertEquals(3, baseDepth.depthAtBase(100));
        assertEquals(2, baseDepth.depthAtBase(115));
        assertEquals(3, baseDepth.depthAtBase(120));

        final Set<Integer> candidateJunctions = Sets.newHashSet(100, 115, 120, 150, 200);
        Map<Integer,Integer> depthMap = baseDepth.createPositionMap(candidateJunctions);

        assertEquals(4, depthMap.size());

        BaseDepth mapDepth = new BaseDepth(baseDepth, depthMap);
        assertEquals(3, mapDepth.depthAtBase(100));
        assertEquals(0, mapDepth.depthAtBase(101));
        assertEquals(2, mapDepth.depthAtBase(115));
        assertEquals(0, mapDepth.depthAtBase(200)); // below the threshold for inclusion
    }
}
