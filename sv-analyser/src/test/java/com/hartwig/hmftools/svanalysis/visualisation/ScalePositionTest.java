package com.hartwig.hmftools.svanalysis.visualisation;

import static org.junit.Assert.assertEquals;

import java.util.Map;
import java.util.Random;

import com.google.common.io.Resources;

import org.junit.Test;

public class ScalePositionTest {

    private static final String LINKS_PATH = Resources.getResource("visualisation/ClusterLinks.tsv").getPath();


    private static final int SCALE_1 = 1;
    private static final int SCALE_10 = 2;
    private static final int SCALE_100 = 5;
    private static final int SCALE_1000 = 10;

//    @Test
//    public void testScaleRegions() throws IOException {
//        final List<GenomeRegion> links = TrackFile.readLinks(LINKS_PATH);
//        final List<GenomeRegion> scaledLinks = ScalePosition.scale(1, links);
//        assertEquals(11, scaledLinks.get(0).start());
//        assertEquals(121, scaledLinks.get(0).end());
//    }

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

}
