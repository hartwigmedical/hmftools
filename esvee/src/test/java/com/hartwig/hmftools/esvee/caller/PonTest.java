package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.utils.PonCombiner.mergeSglRegions;
import static com.hartwig.hmftools.esvee.utils.PonCombiner.mergeSvRegions;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.caller.annotation.PonSglRegion;
import com.hartwig.hmftools.esvee.caller.annotation.PonSvRegion;

import org.junit.Test;

public class PonTest
{
    @Test
    public void testPonSvRegionMerge()
    {
        List<PonSvRegion> combinedRegions = Lists.newArrayList();

        combinedRegions.add(new PonSvRegion(
                new BaseRegion(100, 105), FORWARD,
                new ChrBaseRegion(CHR_2, 100, 110), REVERSE, 10));

        // will be combine with 2
        combinedRegions.add(new PonSvRegion(
                new BaseRegion(98, 103), FORWARD,
                new ChrBaseRegion(CHR_2, 102, 104), REVERSE, 10));

        // not combined, diff orientation
        combinedRegions.add(new PonSvRegion(
                new BaseRegion(102, 107), REVERSE,
                new ChrBaseRegion(CHR_2, 100, 110), REVERSE, 10));

        // not combined, latter breakend doesn't overlap
        combinedRegions.add(new PonSvRegion(
                new BaseRegion(104, 109), FORWARD,
                new ChrBaseRegion(CHR_2, 120, 125), REVERSE, 10));

        // latter start but still overlaps, will be combined
        combinedRegions.add(new PonSvRegion(
                new BaseRegion(105, 106), FORWARD,
                new ChrBaseRegion(CHR_2, 105, 115), REVERSE, 10));

        // new region
        combinedRegions.add(new PonSvRegion(
                new BaseRegion(200, 205), FORWARD,
                new ChrBaseRegion(CHR_2, 100, 105), REVERSE, 10));

        combinedRegions.add(new PonSvRegion(
                new BaseRegion(202, 207), FORWARD,
                new ChrBaseRegion(CHR_2, 95, 100), REVERSE, 10));

        mergeSvRegions(CHR_1, combinedRegions);

        assertEquals(4, combinedRegions.size());

        // check combined regions
        assertEquals(98, combinedRegions.get(0).RegionStart.start());
        assertEquals(106, combinedRegions.get(0).RegionStart.end());
        assertEquals(100, combinedRegions.get(0).RegionEnd.start());
        assertEquals(115, combinedRegions.get(0).RegionEnd.end());

        // unchanged
        assertEquals(102, combinedRegions.get(1).RegionStart.start());
        assertEquals(107, combinedRegions.get(1).RegionStart.end());
        assertEquals(100, combinedRegions.get(1).RegionEnd.start());
        assertEquals(110, combinedRegions.get(1).RegionEnd.end());

        // unchanged
        assertEquals(104, combinedRegions.get(2).RegionStart.start());
        assertEquals(109, combinedRegions.get(2).RegionStart.end());
        assertEquals(120, combinedRegions.get(2).RegionEnd.start());
        assertEquals(125, combinedRegions.get(2).RegionEnd.end());

        // merged
        assertEquals(200, combinedRegions.get(3).RegionStart.start());
        assertEquals(207, combinedRegions.get(3).RegionStart.end());
        assertEquals(95, combinedRegions.get(3).RegionEnd.start());
        assertEquals(105, combinedRegions.get(3).RegionEnd.end());
    }

    @Test
    public void testPonSglRegionMerge()
    {
        List<PonSglRegion> combinedRegions = Lists.newArrayList();

        combinedRegions.add(new PonSglRegion(
                new BaseRegion(100, 105), FORWARD, 10));

        // will be combine with 2
        combinedRegions.add(new PonSglRegion(
                new BaseRegion(98, 103), FORWARD, 10));

        // not combined, diff orientation
        combinedRegions.add(new PonSglRegion(
                new BaseRegion(100, 105), REVERSE, 10));

        // latter start but still overlaps, will be combined
        combinedRegions.add(new PonSglRegion(
                new BaseRegion(105, 106), FORWARD, 10));

        // new region
        combinedRegions.add(new PonSglRegion(
                new BaseRegion(200, 205), REVERSE, 10));

        combinedRegions.add(new PonSglRegion(
                new BaseRegion(203, 210), REVERSE, 10));

        mergeSglRegions(CHR_1, combinedRegions);

        assertEquals(3, combinedRegions.size());

        // check combined regions
        assertEquals(98, combinedRegions.get(0).Region.start());
        assertEquals(106, combinedRegions.get(0).Region.end());

        // unchanged
        assertEquals(100, combinedRegions.get(1).Region.start());
        assertEquals(105, combinedRegions.get(1).Region.end());

        // merged
        assertEquals(200, combinedRegions.get(2).Region.start());
        assertEquals(210, combinedRegions.get(2).Region.end());
    }
}
