package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.RegionUtils.isFullyOverlappedBy;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCentre;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCentreFloat;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionEndingAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionIntersection;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionNegatedIntersection;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionOverlapsOrAdjacent;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionStartingAt;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class RegionUtilsTest
{
    private static final double EPSILON = 1e-9;

    @Test
    public void testRegionNegatedIntersectionEmpty()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        List<BaseRegion> expected = List.of(region);
        List<BaseRegion> actual = regionNegatedIntersection(region, Stream.empty());
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionSingleRegionBeforeTarget()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(new BaseRegion(100, 999));
        List<BaseRegion> expected = List.of(region);
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionMultipleRegionsBeforeTarget()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(
                new BaseRegion(100, 150),
                new BaseRegion(140, 160),
                new BaseRegion(300, 999));
        List<BaseRegion> expected = List.of(region);
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionSingleRegionAfterTarget()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(new BaseRegion(2001, 4000));
        List<BaseRegion> expected = List.of(region);
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionMultipleRegionsAfterTarget()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(
                new BaseRegion(2001, 3000),
                new BaseRegion(2500, 2900),
                new BaseRegion(4000, 5000));
        List<BaseRegion> expected = List.of(region);
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionRegionsOutsideTarget()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(
                new BaseRegion(100, 150),
                new BaseRegion(140, 160),
                new BaseRegion(300, 999),
                new BaseRegion(2001, 3000),
                new BaseRegion(2500, 2900),
                new BaseRegion(4000, 5000));
        List<BaseRegion> expected = List.of(region);
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionRegionsOverlapTargetStart()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(
                new BaseRegion(950, 1050));
        List<BaseRegion> expected = List.of(new BaseRegion(1051, 2000));
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionRegionsOverlapTargetEnd()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(
                new BaseRegion(1050, 2050));
        List<BaseRegion> expected = List.of(new BaseRegion(1000, 1049));
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionSingleRegionWithinTarget()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(
                new BaseRegion(1050, 1080));
        List<BaseRegion> expected = List.of(
                new BaseRegion(1000, 1049),
                new BaseRegion(1081, 2000));
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionOverlappingRegions()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(
                // = start, < end
                new BaseRegion(1010, 1020),
                new BaseRegion(1010, 1015),
                // = start, = end
                new BaseRegion(1030, 1035),
                new BaseRegion(1030, 1035),
                // = start, > end
                new BaseRegion(1050, 1055),
                new BaseRegion(1050, 1060),
                // < start, < end
                new BaseRegion(1110, 1120),
                new BaseRegion(1105, 1115),
                // < start, = end
                new BaseRegion(1130, 1140),
                new BaseRegion(1125, 1140),
                // < start, > end
                new BaseRegion(1150, 1160),
                new BaseRegion(1145, 1165)
        );
        List<BaseRegion> expected = List.of(
                new BaseRegion(1000, 1009),
                new BaseRegion(1021, 1029),
                new BaseRegion(1036, 1049),
                new BaseRegion(1061, 1104),
                new BaseRegion(1121, 1124),
                new BaseRegion(1141, 1144),
                new BaseRegion(1166, 2000));
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionMixed()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(
                new BaseRegion(100, 200),
                new BaseRegion(950, 1050),
                new BaseRegion(1100, 1200),
                new BaseRegion(1070, 1210),
                new BaseRegion(1400, 1700),
                new BaseRegion(1699, 1800));
        List<BaseRegion> expected = List.of(
                new BaseRegion(1051, 1069),
                new BaseRegion(1211, 1399),
                new BaseRegion(1801, 2000)
        );
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionOverlapWholeTarget()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(
                new BaseRegion(1000, 1100),
                new BaseRegion(1101, 1300),
                new BaseRegion(1301, 1400),
                new BaseRegion(1350, 1500),
                new BaseRegion(1500, 1800),
                new BaseRegion(1801, 2000)
        );
        List<BaseRegion> expected = emptyList();
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testRegionNegatedIntersectionTargetWithinRegions()
    {
        BaseRegion region = new BaseRegion(1000, 2000);
        Stream<BaseRegion> regions = Stream.of(new BaseRegion(100, 3000));
        List<BaseRegion> expected = emptyList();
        List<BaseRegion> actual = regionNegatedIntersection(region, regions);
        assertEquals(expected, actual);
    }

    @Test
    public void testIsFullyOverlappedByNoOverlap()
    {
        assertFalse(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.empty()));

        assertFalse(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("2", 10, 20))));
        assertFalse(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("1", 1, 9))));
        assertFalse(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("1", 21, 30))));

        assertFalse(isFullyOverlappedBy(
                new ChrBaseRegion("1", 10, 20),
                Stream.of(
                        new ChrBaseRegion("2", 11, 20),
                        new ChrBaseRegion("1", 25, 30),
                        new ChrBaseRegion("1", 21, 27),
                        new ChrBaseRegion("1", 5, 6),
                        new ChrBaseRegion("10", 11, 20)
                )));
    }

    @Test
    public void testIsFullyOverlappedByPartialOverlap()
    {
        assertFalse(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("1", 1, 10))));
        assertFalse(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("1", 1, 11))));
        assertFalse(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("1", 1, 19))));
        assertFalse(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("1", 20, 30))));
        assertFalse(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("1", 15, 30))));
        assertFalse(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("1", 11, 30))));

        assertFalse(isFullyOverlappedBy(
                new ChrBaseRegion("1", 10, 20),
                Stream.of(
                        new ChrBaseRegion("1", 19, 30),
                        new ChrBaseRegion("1", 13, 15),
                        new ChrBaseRegion("1", 14, 17),
                        new ChrBaseRegion("1", 21, 30),
                        new ChrBaseRegion("1", 5, 11)
                )));
    }

    @Test
    public void testIsFullyOverlappedByFullOverlap()
    {
        assertTrue(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("1", 10, 20))));
        assertTrue(isFullyOverlappedBy(new ChrBaseRegion("1", 10, 20), Stream.of(new ChrBaseRegion("1", 5, 25))));

        assertTrue(isFullyOverlappedBy(
                new ChrBaseRegion("1", 10, 20),
                Stream.of(
                        new ChrBaseRegion("1", 10, 15),
                        new ChrBaseRegion("1", 1, 5),
                        new ChrBaseRegion("1", 19, 23),
                        new ChrBaseRegion("1", 6, 9),
                        new ChrBaseRegion("1", 16, 18),
                        new ChrBaseRegion("1", 5, 10),
                        new ChrBaseRegion("1", 25, 30)
                )));
    }

    @Test
    public void testRegionCentreFloat()
    {
        assertEquals(1.0, regionCentreFloat(new BaseRegion(1, 1)), EPSILON);
        assertEquals(1.5, regionCentreFloat(new BaseRegion(1, 2)), EPSILON);
        assertEquals(5.0, regionCentreFloat(new BaseRegion(1, 9)), EPSILON);
        assertEquals(5.5, regionCentreFloat(new BaseRegion(1, 10)), EPSILON);
        assertEquals(15.0, regionCentreFloat(new BaseRegion(10, 20)), EPSILON);
    }

    @Test
    public void testRegionCentre()
    {
        assertEquals(1, regionCentre(new BaseRegion(1, 1)));
        assertEquals(1, regionCentre(new BaseRegion(1, 2)));
        assertEquals(5, regionCentre(new BaseRegion(1, 9)));
        assertEquals(5, regionCentre(new BaseRegion(1, 10)));
        assertEquals(15, regionCentre(new BaseRegion(10, 20)));
    }

    @Test
    public void testRegionStartingAt()
    {
        assertEquals(new BaseRegion(10, 10), regionStartingAt(10, 1));
        assertEquals(new BaseRegion(1, 10), regionStartingAt(1, 10));
        assertEquals(new BaseRegion(10, 39), regionStartingAt(10, 30));
    }

    @Test
    public void testRegionCenteredAt()
    {
        assertEquals(new BaseRegion(10, 10), regionCenteredAt(10, 1));
        assertEquals(new BaseRegion(10, 11), regionCenteredAt(10, 2));
        assertEquals(new BaseRegion(9, 11), regionCenteredAt(10, 3));
    }

    @Test
    public void testRegionCentreConsistent()
    {
        assertEquals(10, regionCentre(regionCenteredAt(10, 1)));
        assertEquals(10, regionCentre(regionCenteredAt(10, 2)));
        assertEquals(10, regionCentre(regionCenteredAt(10, 3)));
        assertEquals(11, regionCentre(regionCenteredAt(11, 1)));
        assertEquals(11, regionCentre(regionCenteredAt(11, 2)));
        assertEquals(11, regionCentre(regionCenteredAt(11, 3)));
    }

    @Test
    public void testRegionEndingAt()
    {
        assertEquals(new BaseRegion(10, 10), regionEndingAt(10, 1));
        assertEquals(new BaseRegion(1, 10), regionEndingAt(10, 10));
        assertEquals(new BaseRegion(11, 40), regionEndingAt(40, 30));
    }

    @Test
    public void testRegionEndingAtConsistent()
    {
        assertEquals(regionStartingAt(10, 10), regionEndingAt(19, 10));
    }

    @Test
    public void testRegionIntersection()
    {
        assertEquals(
                Optional.empty(),
                regionIntersection(new BaseRegion(10, 20), new BaseRegion(30, 40)));
        assertEquals(
                Optional.empty(),
                regionIntersection(new BaseRegion(10, 20), new BaseRegion(21, 40)));

        assertEquals(
                Optional.of(new BaseRegion(20, 20)),
                regionIntersection(new BaseRegion(10, 20), new BaseRegion(20, 40)));
        assertEquals(
                Optional.of(new BaseRegion(15, 20)),
                regionIntersection(new BaseRegion(10, 20), new BaseRegion(15, 40)));
        assertEquals(
                Optional.of(new BaseRegion(10, 20)),
                regionIntersection(new BaseRegion(10, 20), new BaseRegion(10, 20)));
        assertEquals(
                Optional.of(new BaseRegion(10, 20)),
                regionIntersection(new BaseRegion(10, 20), new BaseRegion(5, 40)));
    }

    @Test
    public void testRegionOverlapsOrAdjacent()
    {
        assertFalse(regionOverlapsOrAdjacent(new BaseRegion(10, 20), new BaseRegion(30, 40)));
        assertFalse(regionOverlapsOrAdjacent(new BaseRegion(30, 40), new BaseRegion(10, 20)));
        assertFalse(regionOverlapsOrAdjacent(new BaseRegion(10, 20), new BaseRegion(22, 30)));

        assertTrue(regionOverlapsOrAdjacent(new BaseRegion(10, 20), new BaseRegion(10, 20)));
        assertTrue(regionOverlapsOrAdjacent(new BaseRegion(10, 20), new BaseRegion(20, 30)));
        assertTrue(regionOverlapsOrAdjacent(new BaseRegion(10, 20), new BaseRegion(21, 30)));
        assertTrue(regionOverlapsOrAdjacent(new BaseRegion(10, 20), new BaseRegion(1, 10)));
        assertTrue(regionOverlapsOrAdjacent(new BaseRegion(10, 20), new BaseRegion(1, 9)));
    }
}
