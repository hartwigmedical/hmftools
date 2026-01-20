package com.hartwig.hmftools.purple.segment;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.junit.Test;

public class SegmentRefinerTest
{
    private final String chr = HumanChromosome._1.shortName();

    @Test
    public void singleRefiningSegmentOneInputIntervals()
    {
        // x is in the middle of a
        PurpleSupportSegment a = pss(101, 200);
        PurpleSupportSegment x = pss(150, 160);
        SegmentRefiner refinement = new SegmentRefiner(List.of(x));
        List<PurpleSupportSegment> result = refinement.refine(List.of(a));
        assertEquals(3, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(x.start() - 1, result.get(0).end());
        assertEquals(x.start(), result.get(1).start());
        assertEquals(x.end(), result.get(1).end());
        assertEquals(x.end() + 1, result.get(2).start());
        assertEquals(a.end(), result.get(2).end());

        // x is at the start of a
        x = pss(90, 110);
        refinement = new SegmentRefiner(List.of(x));
        result = refinement.refine(List.of(a));
        assertEquals(2, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(x.end(), result.get(0).end());
        assertEquals(x.end() + 1, result.get(1).start());
        assertEquals(a.end(), result.get(1).end());

        // x is at the end of a
        x = pss(190, 210);
        refinement = new SegmentRefiner(List.of(x));
        result = refinement.refine(List.of(a));
        assertEquals(2, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(x.start() - 1, result.get(0).end());
        assertEquals(x.start(), result.get(1).start());
        assertEquals(a.end(), result.get(1).end());

        // x does not overlap a
        x = pss(210, 220);
        refinement = new SegmentRefiner(List.of(x));
        result = refinement.refine(List.of(a));
        assertEquals(1, result.size());
        assertEquals(a, result.get(0));
    }

    @Test
    public void singleRefiningSegmentTwoInputIntervals()
    {
        // x straddles a and b
        PurpleSupportSegment a = pss(101, 200);
        PurpleSupportSegment b = pss(201, 300);
        PurpleSupportSegment x = pss(150, 250);
        SegmentRefiner refinement = new SegmentRefiner(List.of(x));
        List<PurpleSupportSegment> result = refinement.refine(List.of(a, b));
        assertEquals(4, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(x.start() - 1, result.get(0).end());
        assertEquals(x.start(), result.get(1).start());
        assertEquals(a.end(), result.get(1).end());
        assertEquals(b.start(), result.get(2).start());
        assertEquals(x.end(), result.get(2).end());
        assertEquals(x.end() + 1, result.get(3).start());
        assertEquals(b.end(), result.get(3).end());

        // start of x is start of b
        x = pss(b.start(), b.end() - 10);
        refinement = new SegmentRefiner(List.of(x));
        result = refinement.refine(List.of(a, b));
        assertEquals(3, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(a.end(), result.get(0).end());
        assertEquals(x.start(), result.get(1).start());
        assertEquals(x.end(), result.get(1).end());
        assertEquals(x.end() + 1, result.get(2).start());
        assertEquals(b.end(), result.get(2).end());

        // end of a is end of x
        x = pss(a.start() + 10, a.end());
        refinement = new SegmentRefiner(List.of(x));
        result = refinement.refine(List.of(a, b));
        assertEquals(3, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(x.start() - 1, result.get(0).end());
        assertEquals(x.start(), result.get(1).start());
        assertEquals(a.end(), result.get(1).end());
        assertEquals(b.start(), result.get(2).start());
        assertEquals(b.end(), result.get(2).end());
    }

    @Test
    public void singleRefiningSegmentThreeInputIntervals()
    {
        PurpleSupportSegment a = pss(101, 200);
        PurpleSupportSegment b = pss(201, 300);
        PurpleSupportSegment c = pss(301, 400);
        PurpleSupportSegment x = pss(150, 350);
        SegmentRefiner refinement = new SegmentRefiner(List.of(x));
        List<PurpleSupportSegment> result = refinement.refine(List.of(a, b, c));
        assertEquals(5, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(x.start() - 1, result.get(0).end());
        assertEquals(x.start(), result.get(1).start());
        assertEquals(a.end(), result.get(1).end());
        assertEquals(b.start(), result.get(2).start());
        assertEquals(b.end(), result.get(2).end());
        assertEquals(c.start(), result.get(3).start());
        assertEquals(x.end(), result.get(3).end());
        assertEquals(x.end() + 1, result.get(4).start());
        assertEquals(c.end(), result.get(4).end());
    }

    @Test
    public void singleRefiningSegmentFourInputIntervals()
    {
        PurpleSupportSegment a = pss(101, 200);
        PurpleSupportSegment b = pss(201, 300);
        PurpleSupportSegment c = pss(301, 400);
        PurpleSupportSegment d = pss(401, 500);
        PurpleSupportSegment x = pss(150, 450);
        SegmentRefiner refinement = new SegmentRefiner(List.of(x));
        List<PurpleSupportSegment> result = refinement.refine(List.of(a, b, c, d));
        assertEquals(6, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(x.start() - 1, result.get(0).end());
        assertEquals(x.start(), result.get(1).start());
        assertEquals(a.end(), result.get(1).end());
        assertEquals(b.start(), result.get(2).start());
        assertEquals(b.end(), result.get(2).end());
        assertEquals(c.start(), result.get(3).start());
        assertEquals(c.end(), result.get(3).end());
        assertEquals(d.start(), result.get(4).start());
        assertEquals(x.end(), result.get(4).end());
        assertEquals(x.end() + 1, result.get(5).start());
        assertEquals(d.end(), result.get(5).end());
    }

    @Test
    public void twoRefiningSegmentsSingleLargeInputInterval()
    {
        PurpleSupportSegment a = pss(101, 200);
        PurpleSupportSegment x = pss(120, 130);
        PurpleSupportSegment y = pss(150, 160);

        SegmentRefiner refinement = new SegmentRefiner(List.of(x, y));
        List<PurpleSupportSegment> result = refinement.refine(List.of(a));
        assertEquals(5, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(x.start() - 1, result.get(0).end());
        assertEquals(x.start(), result.get(1).start());
        assertEquals(x.end(), result.get(1).end());
        assertEquals(x.end() + 1, result.get(2).start());
        assertEquals(y.start() - 1, result.get(2).end());
        assertEquals(y.start(), result.get(3).start());
        assertEquals(y.end(), result.get(3).end());
        assertEquals(y.end() + 1, result.get(4).start());
        assertEquals(a.end(), result.get(4).end());
    }

    @Test
    public void twoRefiningSegmentsTwoIntervals()
    {
        PurpleSupportSegment a = pss(101, 200);
        PurpleSupportSegment b = pss(201, 300);

        // x spans a and b
        PurpleSupportSegment x = pss(180, 220);
        PurpleSupportSegment y = pss(250, 260);

        SegmentRefiner refinement = new SegmentRefiner(List.of(x, y));
        List<PurpleSupportSegment> result = refinement.refine(List.of(a, b));

        assertEquals(6, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(x.start() - 1, result.get(0).end());
        assertEquals(x.start(), result.get(1).start());
        assertEquals(a.end(), result.get(1).end());
        assertEquals(b.start(), result.get(2).start());
        assertEquals(x.end(), result.get(2).end());
        assertEquals(x.end() + 1, result.get(3).start());
        assertEquals(y.start() - 1, result.get(3).end());
        assertEquals(y.start(), result.get(4).start());
        assertEquals(y.end(), result.get(4).end());
        assertEquals(y.end() + 1, result.get(5).start());
        assertEquals(b.end(), result.get(5).end());
    }

    @Test
    public void multipleSegmentsBetweenFilters()
    {
        PurpleSupportSegment a = pss(101, 200);
        PurpleSupportSegment b = pss(201, 300);
        PurpleSupportSegment c = pss(301, 400);
        PurpleSupportSegment d = pss(401, 500);
        PurpleSupportSegment e = pss(501, 600);
        PurpleSupportSegment f = pss(601, 700);
        PurpleSupportSegment g = pss(701, 800);
        PurpleSupportSegment h = pss(801, 900);

        PurpleSupportSegment x = pss(120, 130);
        PurpleSupportSegment y = pss(420, 430);
        PurpleSupportSegment z = pss(780, 820);
        List<GenomeRegion> refiners = List.of(x, y, z);
        List<PurpleSupportSegment> segments = List.of(a, b, c, d, e, f, g, h);

        SegmentRefiner refinement = new SegmentRefiner(refiners);
        List<PurpleSupportSegment> result = refinement.refine(segments);
        assertEquals(14, result.size());
        assertEquals(a.start(), result.get(0).start());
        assertEquals(x.start() - 1, result.get(0).end());
        assertEquals(x.start(), result.get(1).start());
        assertEquals(x.end(), result.get(1).end());
        assertEquals(x.end() + 1, result.get(2).start());
        assertEquals(a.end(), result.get(2).end());
        assertEquals(b, result.get(3));
        assertEquals(c, result.get(4));
        assertEquals(d.start(), result.get(5).start());
        assertEquals(y.start() - 1, result.get(5).end());
        assertEquals(y.start(), result.get(6).start());
        assertEquals(y.end(), result.get(6).end());
        assertEquals(y.end() + 1, result.get(7).start());
        assertEquals(d.end(), result.get(7).end());
        assertEquals(e, result.get(8));
        assertEquals(f, result.get(9));
        assertEquals(g.start(), result.get(10).start());
        assertEquals(z.start() - 1, result.get(10).end());
        assertEquals(z.start(), result.get(11).start());
        assertEquals(g.end(), result.get(11).end());
        assertEquals(h.start(), result.get(12).start());
        assertEquals(z.end(), result.get(12).end());
        assertEquals(z.end() + 1, result.get(13).start());
        assertEquals(h.end(), result.get(13).end());
    }

    @Test
    public void multipleSegmentsBeforeAndAfterFilters()
    {
        PurpleSupportSegment a = pss(101, 200);
        PurpleSupportSegment b = pss(201, 300);
        PurpleSupportSegment c = pss(301, 400);
        PurpleSupportSegment d = pss(401, 500);
        PurpleSupportSegment e = pss(501, 600);
        PurpleSupportSegment f = pss(601, 700);
        PurpleSupportSegment g = pss(701, 800);
        PurpleSupportSegment h = pss(801, 900);

        PurpleSupportSegment x = pss(420, 430);
        PurpleSupportSegment y = pss(440, 450);
        List<GenomeRegion> refiners = List.of(x, y);
        List<PurpleSupportSegment> segments = List.of(a, b, c, d, e, f, g, h);

        SegmentRefiner refinement = new SegmentRefiner(refiners);
        List<PurpleSupportSegment> result = refinement.refine(segments);
        assertEquals(12, result.size());
        assertEquals(a, result.get(0));
        assertEquals(b, result.get(1));
        assertEquals(c, result.get(2));
        assertEquals(d.start(), result.get(3).start());
        assertEquals(x.start() - 1, result.get(3).end());
        assertEquals(x.start(), result.get(4).start());
        assertEquals(x.end(), result.get(4).end());
        assertEquals(x.end() + 1, result.get(5).start());
        assertEquals(y.start() - 1, result.get(5).end());
        assertEquals(y.start(), result.get(6).start());
        assertEquals(y.end(), result.get(6).end());
        assertEquals(y.end() + 1, result.get(7).start());
        assertEquals(d.end(), result.get(7).end());
        assertEquals(e, result.get(8));
        assertEquals(f, result.get(9));
        assertEquals(g, result.get(10));
        assertEquals(h, result.get(11));
    }

    @Test
    public void segmentsTest()
    {
        PurpleSupportSegment x = pss(180, 220);
        PurpleSupportSegment y = pss(250, 160);

        SegmentRefiner refinement = new SegmentRefiner(List.of(x, y));
        List<GenomeRegion> result = refinement.segments();
        assertEquals(2, result.size());
        assertEquals(x, result.get(0));
        assertEquals(y, result.get(1));
    }

    private PurpleSupportSegment pss(int start, int end)
    {
        return new PurpleSupportSegment(chr, start, end, true, SegmentSupport.NONE, false, start, end);
    }
}
