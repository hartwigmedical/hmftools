package com.hartwig.hmftools.svanalysis.visualisation;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class TrackFactoryTest {

    @Test
    public void testTrackStrategies() {
        final String first = "1\t1\t1\t1000";
        final String second = "1\t2\t1\t1000";
        final String third = "1\t1\t1\t1000";

        final List<Track> tracks = TrackFactory.fromString(Lists.newArrayList(first, second, third));

        final List<Track> alwaysIncrement = TrackFactory.alwaysIncrement(tracks);
        assertEquals(3, alwaysIncrement.size());
        assertEquals(1, alwaysIncrement.get(0).track());
        assertEquals(2, alwaysIncrement.get(1).track());
        assertEquals(3, alwaysIncrement.get(2).track());

        final List<Track> incrementOnChromosome = TrackFactory.incrementOnChromosome(tracks);
        assertEquals(3, incrementOnChromosome.size());
        assertEquals(1, incrementOnChromosome.get(0).track());
        assertEquals(1, incrementOnChromosome.get(1).track());
        assertEquals(2, incrementOnChromosome.get(2).track());
    }

}
