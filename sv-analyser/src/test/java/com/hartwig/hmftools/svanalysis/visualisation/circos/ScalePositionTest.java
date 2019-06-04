package com.hartwig.hmftools.svanalysis.visualisation.circos;

import static org.junit.Assert.assertEquals;

import java.util.Map;
import java.util.Random;

import com.google.common.collect.Maps;

import org.junit.Test;

public class ScalePositionTest {

    private static final int SCALE_1 = 1;
    private static final int SCALE_10 = 2;
    private static final int SCALE_100 = 9;
    private static final int SCALE_1000 = 28;

    @Test
    public void testFirstPositionIsAtStart() {
        int start = new Random().nextInt(1000);
        long firstPosition = new Random().nextInt(1000000) + 1;

        assertEquals(start, (int) ScalePosition.positionMap(start, firstPosition).get(firstPosition));
    }

    @Test
    public void testScalePosition() {
        int start = new Random().nextInt(1000);
        long firstPosition = new Random().nextInt(1000000) + 1;
        final Map<Long, Integer> map =
                ScalePosition.positionMap(start, firstPosition + 1001, firstPosition, firstPosition + 1, firstPosition);
        assertEquals(3, map.size());

        assertEquals(start, (int) map.get(firstPosition));
        assertEquals(map.get(firstPosition) + SCALE_1, (int) map.get(firstPosition + 1));
        assertEquals(map.get(firstPosition + 1) + SCALE_1000, (int) map.get(firstPosition + 1001));
    }

    @Test
    public void testLogDistance() {
        assertEquals(SCALE_1, ScalePosition.logDistance(1));
        assertEquals(SCALE_10, ScalePosition.logDistance(10));
        assertEquals(SCALE_100, ScalePosition.logDistance(100));
        assertEquals(SCALE_1000, ScalePosition.logDistance(1000));
    }

    @Test
    public void testInterpolate() {
        final Map<Long, Integer> map = Maps.newHashMap();
        map.put(1000L, 10);

        assertEquals(10, ScalePosition.interpolate(900L, map));
        assertEquals(10, ScalePosition.interpolate(1000L, map));
        assertEquals(10, ScalePosition.interpolate(1100L, map));

        map.put(2000L, 20);
        assertEquals(10, ScalePosition.interpolate(900L, map));
        assertEquals(10, ScalePosition.interpolate(1000L, map));
        assertEquals(11, ScalePosition.interpolate(1100L, map));
        assertEquals(15, ScalePosition.interpolate(1500L, map));
        assertEquals(20, ScalePosition.interpolate(2000L, map));
        assertEquals(20, ScalePosition.interpolate(2100L, map));

        map.put(3000L, 40);
        assertEquals(20, ScalePosition.interpolate(2000L, map));
        assertEquals(22, ScalePosition.interpolate(2100L, map));
        assertEquals(30, ScalePosition.interpolate(2500L, map));
        assertEquals(40, ScalePosition.interpolate(3000L, map));
    }

}
