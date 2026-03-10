package com.hartwig.hmftools.amber.purity;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class RegionsFilterTest
{
    List<GenomePositionImpl> Positions = new ArrayList<>();

    @Test
    public void multipleRegionsTest()
    {
        List<ChrBaseRegion> regions = new ArrayList<>();
        regions.add(region("3", 40_001, 50_000));
        regions.add(region("2", 40_001, 50_000));
        regions.add(region("1", 40_001, 50_000));

        regions.add(region("2", 10_001, 20_000));
        regions.add(region("3", 10_001, 20_000));
        regions.add(region("1", 10_001, 20_000));

        regions.add(region("2", 1001, 2000));
        regions.add(region("3", 1001, 2000));
        regions.add(region("1", 1001, 2000));

        RegionsFilter filter = new RegionsFilter(regions);

        addPosition("1", 100); // in
        addPosition("1", 200); // in
        addPosition("1", 900); // in
        addPosition("1", 1000); // in
        addPosition("1", 1001);
        addPosition("1", 1200);
        addPosition("1", 1300);
        addPosition("1", 1900);
        addPosition("1", 2000);
        addPosition("1", 2001); // in
        addPosition("1", 20_001); // in
        addPosition("1", 30_001); // in
        addPosition("1", 40_001);
        addPosition("1", 50_001); // in
        addPosition("1", 60_001); // in

        addPosition("3", 500); // in
        addPosition("3", 1500);
        addPosition("3", 2500); // in
        addPosition("3", 9500); // in
        addPosition("3", 10_500);
        addPosition("3", 15_500);
        addPosition("3", 20_500); // in
        addPosition("3", 30_500); // in
        addPosition("3", 40_500);
        addPosition("3", 50_500); // in

        addPosition("5", 2500); // in

        List<GenomePositionImpl> filtered = filter.filter(Positions);
        assertEquals(16, filtered.size());
        assertEquals(Positions.get(0), filtered.get(0));
        assertEquals(Positions.get(1), filtered.get(1));
        assertEquals(Positions.get(3), filtered.get(3));
        assertEquals(Positions.get(9), filtered.get(4));
        assertEquals(Positions.get(10), filtered.get(5));
        assertEquals(Positions.get(11), filtered.get(6));
        assertEquals(Positions.get(13), filtered.get(7));
        assertEquals(Positions.get(14), filtered.get(8));
        assertEquals(Positions.get(15), filtered.get(9));
        assertEquals(Positions.get(17), filtered.get(10));
        assertEquals(Positions.get(18), filtered.get(11));
        assertEquals(Positions.get(21), filtered.get(12));
        assertEquals(Positions.get(22), filtered.get(13));
        assertEquals(Positions.get(24), filtered.get(14));
        assertEquals(Positions.get(25), filtered.get(15));
    }

    private ChrBaseRegion region(String chromosome, int start, int end)
    {
        return new ChrBaseRegion(chromosome, start, end);
    }

    private void addPosition(String chromosome, int position)
    {
        Positions.add(new GenomePositionImpl(chromosome, position));
    }
}
