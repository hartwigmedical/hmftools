package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.results.TranscriptResult.calcEffectiveLength;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class MiscTest
{
    @Test
    public void testEffectiveLength()
    {
        final List<int[]> fragmentLengthData = Lists.newArrayList();

        fragmentLengthData.add(new int[]{100, 1});
        fragmentLengthData.add(new int[]{200, 2});
        fragmentLengthData.add(new int[]{300, 1});

        assertEquals(800, calcEffectiveLength(1000, fragmentLengthData), 0.001);

        assertEquals(83.33, calcEffectiveLength(250, fragmentLengthData), 0.1);
    }
}
