package com.hartwig.hmftools.linx.visualiser.circos;

import static org.junit.Assert.assertEquals;

import java.util.Map;
import java.util.Random;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.junit.Test;

public class ScaleContigTest
{

    private static final int SCALE_1 = 10;
    private static final int SCALE_10 = 11;
    private static final int SCALE_100 = 18;
    private static final int SCALE_1000 = 37;

    @Test
    public void testFirstPositionIsAtStart()
    {
        ScaleContig victim = new ScaleContig("1", Lists.newArrayList(200L, 10L));
        assertEquals(1, victim.scale(10L));
    }


    @Test
    public void testScalePosition()
    {
        long firstPosition = new Random().nextInt(1000000) + 1;
        ScaleContig victim =
                new ScaleContig("1", Lists.newArrayList(firstPosition + 1001, firstPosition, firstPosition + 1, firstPosition));

        assertEquals(1, victim.scale(firstPosition));
        assertEquals(1 + SCALE_1, victim.scale(firstPosition + 1));
        assertEquals(1 + SCALE_1 + SCALE_1000, victim.scale(firstPosition + 1001));
    }

    @Test
    public void testLogDistance()
    {
        assertEquals(SCALE_1, ScaleContig.logDistance(1));
        assertEquals(SCALE_10, ScaleContig.logDistance(10));
        assertEquals(SCALE_100, ScaleContig.logDistance(100));
        assertEquals(SCALE_1000, ScaleContig.logDistance(1000));
    }

    @Test
    public void testInterpolate()
    {
        final Map<Long, Integer> map = Maps.newHashMap();
        map.put(1000L, 10);
        final ScaleContig victim = new ScaleContig("Dummy", map);

        assertEquals(10, victim.interpolate(900L));
        assertEquals(10, victim.interpolate(1000L));
        assertEquals(10, victim.interpolate(1100L));

        map.put(2000L, 20);
        assertEquals(10, victim.interpolate(900L));
        assertEquals(10, victim.interpolate(1000L));
        assertEquals(11, victim.interpolate(1100L));
        assertEquals(15, victim.interpolate(1500L));
        assertEquals(20, victim.interpolate(2000L));
        assertEquals(20, victim.interpolate(2100L));

        map.put(3000L, 40);
        assertEquals(20, victim.interpolate(2000L));
        assertEquals(22, victim.interpolate(2100L));
        assertEquals(30, victim.interpolate(2500L));
        assertEquals(40, victim.interpolate(3000L));
    }

}
