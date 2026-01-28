package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.purple.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.SegmentSupport.EXCL;
import static com.hartwig.hmftools.common.purple.SegmentSupport.INV;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionImpl;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.junit.Test;

public class PurpleSupportSegmentTest
{
    private final String chr = HumanChromosome._1.shortName();

    @Test
    public void supportWhenSplit()
    {
        PurpleSupportSegment a = new PurpleSupportSegment(chr, 100, 200, true, CENTROMERE, false, 100, 200);
        List<PurpleSupportSegment> result = a.split(150);
        assertEquals(2, result.size());
        assertEquals(CENTROMERE, result.get(0).Support);
        assertEquals(EXCL, result.get(1).Support);

        a = new PurpleSupportSegment(chr, 100, 200, true, INV, false, 100, 200);
        result = a.split(150);
        assertEquals(2, result.size());
        assertEquals(INV, result.get(0).Support);
        assertEquals(EXCL, result.get(1).Support);

        a = new PurpleSupportSegment(chr, 100, 200, true, EXCL, false, 100, 200);
        result = a.split(150);
        assertEquals(2, result.size());
        assertEquals(EXCL, result.get(0).Support);
        assertEquals(EXCL, result.get(1).Support);
    }

    @Test
    public void supportWhenSplitBySegment()
    {
        PurpleSupportSegment a = new PurpleSupportSegment(chr, 100, 200, true, CENTROMERE, false, 100, 200);
        PurpleSupportSegment b = pss(150, 180);
        List<PurpleSupportSegment> result = a.splitBy(b);
        assertEquals(3, result.size());

        assertEquals(CENTROMERE, result.get(0).Support);
        assertEquals(EXCL, result.get(1).Support);
        assertEquals(EXCL, result.get(2).Support);
    }

    @Test
    public void splitAtPositionBeforeStart()
    {
        PurpleSupportSegment segment = pss(100, 200);
        List<PurpleSupportSegment> result = segment.split(50);
        assertEquals(1, result.size());
        assertEquals(segment, result.get(0));
    }

    @Test
    public void splitAtStart()
    {
        PurpleSupportSegment segment = pss(100, 200);
        List<PurpleSupportSegment> result = segment.split(100);
        assertEquals(1, result.size());
        assertEquals(chr, result.get(0).chromosome());
        assertEquals(100, result.get(0).start());
        assertEquals(200, result.get(0).end());
    }

    @Test
    public void splitWithin()
    {
        PurpleSupportSegment segment = pss(100, 200);
        List<PurpleSupportSegment> result = segment.split(150);
        assertEquals(2, result.size());
        assertEquals(chr, result.get(0).chromosome());
        assertEquals(100, result.get(0).start());
        assertEquals(149, result.get(0).end());
        assertEquals(chr, result.get(1).chromosome());
        assertEquals(150, result.get(1).start());
        assertEquals(200, result.get(1).end());
    }

    @Test
    public void splitAtEnd()
    {
        PurpleSupportSegment segment = pss(100, 200);
        List<PurpleSupportSegment> result = segment.split(200);
        assertEquals(2, result.size());
        assertEquals(chr, result.get(0).chromosome());
        assertEquals(100, result.get(0).start());
        assertEquals(199, result.get(0).end());
        assertEquals(chr, result.get(1).chromosome());
        assertEquals(200, result.get(1).start());
        assertEquals(200, result.get(1).end());
    }

    @Test
    public void splitAtPositionAfterEnd()
    {
        PurpleSupportSegment segment = pss(100, 200);
        List<PurpleSupportSegment> result = segment.split(250);
        assertEquals(1, result.size());
        assertEquals(segment, result.get(0));
    }

    @Test
    public void minAndMaxStartWhenSplitBeforeMaxStart()
    {
        PurpleSupportSegment segment = new PurpleSupportSegment(chr, 100, 200, true, SegmentSupport.NONE, false, 80, 120);
        // Split before the max start
        List<PurpleSupportSegment> result = segment.split(110);
        assertEquals(2, result.size());
        assertEquals(100, result.get(0).start());
        assertEquals(109, result.get(0).end());
        assertEquals(80, result.get(0).minStart());
        assertEquals(109, result.get(0).maxStart());
        assertEquals(chr, result.get(1).chromosome());
        assertEquals(110, result.get(1).start());
        assertEquals(200, result.get(1).end());
        assertEquals(110, result.get(1).minStart());
        assertEquals(110, result.get(1).maxStart());
    }

    @Test
    public void minAndMaxStartWhenSplitAfterMaxStart()
    {
        PurpleSupportSegment segment = new PurpleSupportSegment(chr, 100, 200, true, SegmentSupport.NONE, false, 80, 120);
        // Split before the max start
        List<PurpleSupportSegment> result = segment.split(130);
        assertEquals(2, result.size());
        assertEquals(100, result.get(0).start());
        assertEquals(129, result.get(0).end());
        assertEquals(80, result.get(0).minStart());
        assertEquals(120, result.get(0).maxStart());
        assertEquals(chr, result.get(1).chromosome());
        assertEquals(130, result.get(1).start());
        assertEquals(200, result.get(1).end());
        assertEquals(130, result.get(1).minStart());
        assertEquals(130, result.get(1).maxStart());
    }

    @Test
    public void splitByWithNoIntersection()
    {
        PurpleSupportSegment a = pss(100, 200);
        PurpleSupportSegment b = pss(250, 350);
        List<PurpleSupportSegment> result = a.splitBy(b);
        assertEquals(1, result.size());
        assertEquals(a, result.get(0));
    }

    @Test
    public void splitByWithIntervalOnDifferentChromosome()
    {
        PurpleSupportSegment a = pss(100, 200);
        PurpleSupportSegment b = new PurpleSupportSegment(chr + "2", 150, 250, true, SegmentSupport.NONE, false, 100, 100);
        List<PurpleSupportSegment> result = a.splitBy(b);
        assertEquals(1, result.size());
    }

    @Test
    public void splitByIdenticalSegment()
    {
        PurpleSupportSegment a = pss(100, 200);
        PurpleSupportSegment b = pss(100, 200);
        List<PurpleSupportSegment> result = a.splitBy(b);
        assertEquals(1, result.size());
        assertEquals(a, result.get(0));
    }

    @Test
    public void splitByWithIntervalThatIntersectsAtEnd()
    {
        PurpleSupportSegment a = pss(100, 200);
        PurpleSupportSegment b = pss(150, 250);
        List<PurpleSupportSegment> result = a.splitBy(b);
        assertEquals(2, result.size());
        assertEquals(chr, result.get(0).chromosome());

        assertEquals(a.start(), result.get(0).start());
        assertEquals(b.start() - 1, result.get(0).end());

        assertEquals(b.start(), result.get(1).start());
        assertEquals(a.end(), result.get(1).end());
    }

    @Test
    public void fragmentMinAndMaxStartWhenIntersectionIsAtEnd()
    {
        // First with a.maxStart < b.start
        PurpleSupportSegment a = new PurpleSupportSegment(chr, 100, 200, true, SegmentSupport.NONE, false, 80, 120);
        GenomeRegion b = pss(150, 250);
        List<PurpleSupportSegment> result = a.splitBy(b);
        assertEquals(2, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(a.minStart(), result.get(0).minStart());
        assertEquals(a.maxStart(), result.get(0).maxStart());

        assertEquals(b.start(), result.get(1).start());
        assertEquals(b.start(), result.get(1).minStart());
        assertEquals(b.start(), result.get(1).maxStart());

        // Now with a.maxStart > b.start
        a = new PurpleSupportSegment(chr, 100, 200, true, SegmentSupport.NONE, false, 80, 180);
        result = a.splitBy(b);
        assertEquals(2, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(b.start() - 1, result.get(0).end());
        assertEquals(a.minStart(), result.get(0).minStart());
        assertEquals(result.get(0).end(), result.get(0).maxStart());

        assertEquals(b.start(), result.get(1).start());
        assertEquals(b.start(), result.get(1).minStart());
        assertEquals(b.start(), result.get(1).maxStart());
    }

    @Test
    public void splitBySubInterval()
    {
        PurpleSupportSegment a = pss(100, 200);
        PurpleSupportSegment b = pss(150, 180);
        List<PurpleSupportSegment> result = a.splitBy(b);
        assertEquals(3, result.size());

        assertEquals(a.start(), result.get(0).start());
        assertEquals(b.start() - 1, result.get(0).end());

        assertEquals(b.start(), result.get(1).start());
        assertEquals(b.end(), result.get(1).end());

        assertEquals(b.end() + 1, result.get(2).start());
        assertEquals(a.end(), result.get(2).end());
    }

    @Test
    public void splitByIntervalThatIntersectsAtStart()
    {
        PurpleSupportSegment a = pss(100, 200);
        PurpleSupportSegment b = pss(50, 180);
        List<PurpleSupportSegment> result = a.splitBy(b);
        assertEquals(2, result.size());

        assertEquals(a.start(), result.get(0).start());
        assertEquals(b.end(), result.get(0).end());
        assertEquals(b.end() + 1, result.get(1).start());
        assertEquals(a.end(), result.get(1).end());
    }

    @Test
    public void minAndMaxStartWhenIntersectionIsAtStart()
    {
        // First with a.maxStart < b.end
        PurpleSupportSegment a = new PurpleSupportSegment(chr, 100, 200, true, SegmentSupport.NONE, false, 80, 120);
        GenomeRegion b = new GenomeRegionImpl(chr, 50, 180);
        List<PurpleSupportSegment> result = a.splitBy(b);
        assertEquals(2, result.size());

        assertEquals(a.start(), result.get(0).start());
        assertEquals(a.minStart(), result.get(0).minStart());
        assertEquals(a.maxStart(), result.get(0).maxStart());
        assertEquals(b.end() + 1, result.get(1).start());
        assertEquals(b.end() + 1, result.get(1).minStart());
        assertEquals(b.end() + 1, result.get(1).maxStart());

        // Now with a.maxStart > b.end
        a = new PurpleSupportSegment(chr, 100, 200, true, SegmentSupport.NONE, false, 80, 200);
        result = a.splitBy(b);
        assertEquals(a.start(), result.get(0).start());
        assertEquals(b.end(), result.get(0).end());
        assertEquals(a.minStart(), result.get(0).minStart());
        assertEquals(b.end(), result.get(0).maxStart());
        assertEquals(b.end() + 1, result.get(1).start());
        assertEquals(b.end() + 1, result.get(1).minStart());
        assertEquals(b.end() + 1, result.get(1).maxStart());
    }

    @Test
    public void differenceWithContainingInterval()
    {
        PurpleSupportSegment a = pss(100, 200);
        PurpleSupportSegment b = pss(50, 250);
        List<PurpleSupportSegment> result = a.splitBy(b);
        assertEquals(1, result.size());
        assertEquals(a, result.get(0));
    }

    private PurpleSupportSegment pss(int start, int end)
    {
        return new PurpleSupportSegment(chr, start, end, true, SegmentSupport.NONE, false, start, end);
    }
}
