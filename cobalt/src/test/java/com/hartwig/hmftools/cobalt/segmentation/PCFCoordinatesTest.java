package com.hartwig.hmftools.cobalt.segmentation;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.segmentation.PiecewiseConstantFit;

import org.junit.Test;

public class PCFCoordinatesTest
{
    private String chromosome = "1";

    @Test
    public void intervalsTest()
    {
        chromosome = "1";
        PiecewiseConstantFit pcf =
                new PiecewiseConstantFit(new int[] { 41, 30, 30 }, new int[] { 0, 41, 71 }, new double[] { 0.0, 0.0, 0.0 });
        PCFCoordinates coordinates = new PCFCoordinates(pcf, 1000, 0, chromosome);
        List<ChrBaseRegion> intervals = coordinates.intervals();
        assertEquals(3, intervals.size());
        assertEquals(cbr(0, 40999), intervals.get(0));
        assertEquals(cbr(41000, 70999), intervals.get(1));
        assertEquals(cbr(71000, 100999), intervals.get(2));
    }

    @Test
    public void intervalsWithGapsTest()
    {
        chromosome = "chr2";
        PiecewiseConstantFit pcf =
                new PiecewiseConstantFit(new int[] { 4, 20, 12 }, new int[] { 100, 141, 271 }, new double[] { 0.0, 0.0, 0.0 });
        PCFCoordinates coordinates = new PCFCoordinates(pcf, 1000, 1, chromosome);
        List<ChrBaseRegion> intervals = coordinates.intervals();
        assertEquals(3, intervals.size());
        assertEquals(cbr(100_001, 104000), intervals.get(0));
        assertEquals(cbr(14_1001, 161000), intervals.get(1));
        assertEquals(cbr(271_001, 283000), intervals.get(2));
    }

    private ChrBaseRegion cbr(int start, int end)
    {
        return new ChrBaseRegion(chromosome, start, end);
    }
}
