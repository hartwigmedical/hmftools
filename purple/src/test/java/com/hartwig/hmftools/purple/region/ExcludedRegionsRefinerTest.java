package com.hartwig.hmftools.purple.region;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.purple.segment.PurpleSupportSegment;
import com.hartwig.hmftools.purple.segment.SegmentRefiner;

import org.junit.Before;
import org.junit.Test;

public class ExcludedRegionsRefinerTest
{
    private final String chr1 = "1";
    private final String chr2 = "2";
    private final String chr6 = "6";
    private final String chr7 = "7";
    private final String chr14 = "14";
    Set<ChrBaseRegion> regions = new HashSet<>();
    ExcludedRegionsRefiner refiner = new ExcludedRegionsRefiner(regions);

    @Before
    public void setUp()
    {
        regions = new HashSet<>();
        regions.add(region(chr6, 2000, 2200));
        regions.add(region(chr7, 100, 200));
        regions.add(region(chr6, 1000, 1200));
        regions.add(region(chr14, 2000, 2200));
        regions.add(region(chr6, 100, 200));
        regions.add(region(chr14, 1000, 1200));
        refiner = new ExcludedRegionsRefiner(regions);
    }

    @Test
    public void hasExcludedRegionsTest()
    {
        assert refiner.hasExcludedRegions(chr6);
        assert refiner.hasExcludedRegions(chr7);
        assert refiner.hasExcludedRegions(chr14);
        assert !refiner.hasExcludedRegions(chr1);
        assert !refiner.hasExcludedRegions(chr2);
    }

    @Test
    public void refinerForChromosomeTest()
    {
        SegmentRefiner chr6Refiner = refiner.refiner(chr6);
        List<GenomeRegion> regions = chr6Refiner.segments();
        assertEquals(3, regions.size());
        assertEquals(region(chr6, 100, 200).genomeRegion(), regions.get(0));
        assertEquals(region(chr6, 1000, 1200).genomeRegion(), regions.get(1));
        assertEquals(region(chr6, 2000, 2200).genomeRegion(), regions.get(2));
    }

    @Test
    public void refineTest()
    {
        List<PurpleSupportSegment> inputs = new ArrayList<>();
        inputs.add(segment(chr1, 1001, 2000));
        inputs.add(segment(chr1, 2001, 3000));
        inputs.add(segment(chr1, 3001, 4000));
        inputs.add(segment(chr1, 4001, 5000));
        inputs.add(segment(chr2, 1001, 2000));
        inputs.add(segment(chr2, 2001, 3000));
        inputs.add(segment(chr2, 3001, 4000));
        inputs.add(segment(chr2, 4001, 5000));
        inputs.add(segment(chr6, 1001, 2000));  // becomes [1001, 1200],[1201,1999],[2000,2000]
        inputs.add(segment(chr6, 2001, 3000));  // becomes [2001, 2200],[2201,3000]
        inputs.add(segment(chr6, 3001, 4000));
        inputs.add(segment(chr6, 4001, 5000));
        inputs.add(segment(chr7, 1001, 2000));
        inputs.add(segment(chr7, 2001, 3000));
        inputs.add(segment(chr7, 3001, 4000));
        inputs.add(segment(chr7, 4001, 5000));

        List<PurpleSupportSegment> outputs = refiner.refine(inputs);
        assertEquals(19, outputs.size());
        assertEquals(inputs.get(0), outputs.get(0));
        assertEquals(inputs.get(1), outputs.get(1));
        assertEquals(inputs.get(7), outputs.get(7));
        checkEqual(segment(chr6, 1001, 1200), outputs.get(8));
        checkEqual(segment(chr6, 1201, 1999), outputs.get(9));
        checkEqual(segment(chr6, 2000, 2000), outputs.get(10));
        checkEqual(segment(chr6, 2001, 2200), outputs.get(11));
        checkEqual(segment(chr6, 2201, 3000), outputs.get(12));
        assertEquals(inputs.get(10), outputs.get(13));
        assertEquals(inputs.get(15), outputs.get(18));
    }

    private void checkEqual(PurpleSupportSegment expected, PurpleSupportSegment actual)
    {
        assertEquals(expected.chromosome(), actual.chromosome());
        assertEquals(expected.start(), actual.start());
        assertEquals(expected.end(), actual.end());
    }

    private ChrBaseRegion region(final String chromosome, final int start, final int end)
    {
        return new ChrBaseRegion(chromosome, start, end);
    }

    private PurpleSupportSegment segment(final String chromosome, final int start, final int end)
    {
        return new PurpleSupportSegment(chromosome, start, end, false, SegmentSupport.NONE, false, start, end);
    }
}
