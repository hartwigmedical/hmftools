package com.hartwig.hmftools.svanalysis.visualisation.data;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class SegmentsTest {

    @Test
    public void testTrackStrategies() {
        final String first = "sample,1,1,1,1,1000,1";
        final String second = "sample,1,1,2,1,1000,1";
        final String third = "sample,1,1,1,1,1000,1";

        final List<Segment> segments = Segments.fromString(Lists.newArrayList(first, second, third));

        final List<Segment> alwaysIncrement = Segments.alwaysIncrement(segments);
        assertEquals(3, alwaysIncrement.size());
        assertEquals(1, alwaysIncrement.get(0).track());
        assertEquals(2, alwaysIncrement.get(1).track());
        assertEquals(3, alwaysIncrement.get(2).track());

        final List<Segment> incrementOnChromosome = Segments.incrementOnChromosome(segments);
        assertEquals(3, incrementOnChromosome.size());
        assertEquals(1, incrementOnChromosome.get(0).track());
        assertEquals(1, incrementOnChromosome.get(1).track());
        assertEquals(2, incrementOnChromosome.get(2).track());
    }

}
