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

        /*
        fragmentLengthData.clear();
        fragmentLengthData.add(new int[]{50,31679});
        fragmentLengthData.add(new int[]{100, 266724});
        fragmentLengthData.add(new int[]{150,273855});
        fragmentLengthData.add(new int[]{200,177957});
        fragmentLengthData.add(new int[]{250,107339});
        fragmentLengthData.add(new int[]{300,109962});
        fragmentLengthData.add(new int[]{550,12691});

        double length = calcEffectiveLength(8415, fragmentLengthData);
        // exonicBases(8415) neg calculated effective length(-523.5854130811146) from dist()
        */
    }
}
