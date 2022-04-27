package com.hartwig.hmftools.linx.visualiser.data;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.visualiser.circos.SegmentTerminal;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;

import org.junit.Test;

public class SegmentsTest
{
    @Test
    public void testTrackStrategies()
    {
        final List<VisSegment> segments = Lists.newArrayList(
                createSegment("1"),
                createSegment("2"), createSegment("1"));

        final List<VisSegment> incrementOnChromosome = VisSegments.incrementOnChromosome(segments, Collections.emptyList(), false);
        assertEquals(3, incrementOnChromosome.size());
        assertEquals(1, incrementOnChromosome.get(0).Track);
        assertEquals(1, incrementOnChromosome.get(1).Track);
        assertEquals(2, incrementOnChromosome.get(2).Track);
    }

    private static VisSegment createSegment(String chromosome)
    {
        VisSegment segment = new VisSegment(
                "sampleId", 1,1, chromosome, "1", "1000", 0, false);

        segment.setTerminalStart(SegmentTerminal.NONE);
        segment.setTerminalEnd(SegmentTerminal.NONE);
        segment.Track = 1;
        return segment;
    }

}
