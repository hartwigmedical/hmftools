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
        ScaleContig victim = new ScaleContig("1", Lists.newArrayList(200, 10));
        assertEquals(1, victim.scale(10));
    }


    @Test
    public void testScalePosition()
    {
        int firstPosition = new Random().nextInt(1000000) + 1;
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
        final Map<Integer,Integer> map = Maps.newHashMap();
        map.put(1000, 10);
        final ScaleContig victim = new ScaleContig("Dummy", map);

        assertEquals(10, victim.interpolate(900));
        assertEquals(10, victim.interpolate(1000));
        assertEquals(10, victim.interpolate(1100));

        map.put(2000, 20);
        assertEquals(10, victim.interpolate(900));
        assertEquals(10, victim.interpolate(1000));
        assertEquals(11, victim.interpolate(1100));
        assertEquals(15, victim.interpolate(1500));
        assertEquals(20, victim.interpolate(2000));
        assertEquals(20, victim.interpolate(2100));

        map.put(3000, 40);
        assertEquals(20, victim.interpolate(2000));
        assertEquals(22, victim.interpolate(2100));
        assertEquals(30, victim.interpolate(2500));
        assertEquals(40, victim.interpolate(3000));
    }

}
