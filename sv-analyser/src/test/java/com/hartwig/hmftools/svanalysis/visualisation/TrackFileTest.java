package com.hartwig.hmftools.svanalysis.visualisation;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class TrackFileTest {

    @Test
    public void testTrackStrategies() {
        final String first = "1\t1\t1\t1000";
        final String second = "1\t2\t1\t1000";
        final String third = "1\t1\t1\t1000";

        final List<String> lines = Lists.newArrayList(first, second, third);

        final List<Track> alwaysBump = TrackFile.alwaysBumpTrack(lines);
        assertEquals(3, alwaysBump.size());
        assertEquals(1, alwaysBump.get(0).track());
        assertEquals(2, alwaysBump.get(1).track());
        assertEquals(3, alwaysBump.get(2).track());

        final List<Track> maintainTrackBetweenChromosomes = TrackFile.maintainTrackBetweenChromosomes(lines);
        assertEquals(3, maintainTrackBetweenChromosomes.size());
        assertEquals(1, maintainTrackBetweenChromosomes.get(0).track());
        assertEquals(1, maintainTrackBetweenChromosomes.get(1).track());
        assertEquals(2, maintainTrackBetweenChromosomes.get(2).track());
    }

}
