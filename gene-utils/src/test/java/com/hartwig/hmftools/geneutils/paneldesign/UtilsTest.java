package com.hartwig.hmftools.geneutils.paneldesign;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BaseRegion;

import org.junit.Test;

public class UtilsTest
{
    @Test
    public void testComputeUncoveredRegionsEmptyCovered()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        List<BaseRegion> expected = List.of(target);
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, Stream.empty());
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsSingleCoveredBeforeTarget()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(new BaseRegion(100, 999));
        List<BaseRegion> expected = List.of(target);
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsMultipleCoveredBeforeTarget()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(
                new BaseRegion(100, 150),
                new BaseRegion(140, 160),
                new BaseRegion(300, 999));
        List<BaseRegion> expected = List.of(target);
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsSingleCoveredAfterTarget()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(new BaseRegion(2001, 4000));
        List<BaseRegion> expected = List.of(target);
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsMultipleCoveredAfterTarget()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(
                new BaseRegion(2001, 3000),
                new BaseRegion(2500, 2900),
                new BaseRegion(4000, 5000));
        List<BaseRegion> expected = List.of(target);
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsCoveredOutsideTarget()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(
                new BaseRegion(100, 150),
                new BaseRegion(140, 160),
                new BaseRegion(300, 999),
                new BaseRegion(2001, 3000),
                new BaseRegion(2500, 2900),
                new BaseRegion(4000, 5000));
        List<BaseRegion> expected = List.of(target);
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsCoveredOverlapTargetStart()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(
                new BaseRegion(950, 1050));
        List<BaseRegion> expected = List.of(new BaseRegion(1051, 2000));
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsCoveredOverlapTargetEnd()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(
                new BaseRegion(1050, 2050));
        List<BaseRegion> expected = List.of(new BaseRegion(1000, 1049));
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsSingleCoveredWithinTarget()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(
                new BaseRegion(1050, 1080));
        List<BaseRegion> expected = List.of(
                new BaseRegion(1000, 1049),
                new BaseRegion(1081, 2000));
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsOverlappingCovered()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(
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
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsMixed()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(
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
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsCoveredWholeTarget()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(
                new BaseRegion(1000, 1100),
                new BaseRegion(1101, 1300),
                new BaseRegion(1301, 1400),
                new BaseRegion(1350, 1500),
                new BaseRegion(1500, 1800),
                new BaseRegion(1801, 2000)
        );
        List<BaseRegion> expected = Collections.emptyList();
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testComputeUncoveredRegionsTargetWithinCovered()
    {
        BaseRegion target = new BaseRegion(1000, 2000);
        Stream<BaseRegion> covered = Stream.of(new BaseRegion(100, 3000));
        List<BaseRegion> expected = Collections.emptyList();
        List<BaseRegion> actual = Utils.computeUncoveredRegions(target, covered);
        assertEquals(expected, actual);
    }

    @Test
    public void testGetBestScoringElementEmptyStream()
    {
        assertEquals(Optional.empty(), Utils.<Double>getBestScoringElement(Stream.empty(), s -> s, s -> s == 0, true));

        assertEquals(Optional.empty(), Utils.<Double>getBestScoringElement(Stream.empty(), s -> s, s -> s == 0, false));
    }

    @Test
    public void testGetBestScoringElementOptimal()
    {
        assertEquals(Optional.of(3), Utils.getBestScoringElement(Stream.of(1, 2, 3, 4, 5), s -> s, s -> s == 3, true));
        assertEquals(Optional.of(3), Utils.getBestScoringElement(Stream.of(3, 2, 1, 4, 5), s -> s, s -> s == 3, true));

        assertEquals(Optional.of(3), Utils.getBestScoringElement(Stream.of(5, 4, 3, 2, 1), s -> s, s -> s == 3, false));
        assertEquals(Optional.of(3), Utils.getBestScoringElement(Stream.of(3, 4, 5, 2, 1), s -> s, s -> s == 3, false));
    }

    @Test
    public void testGetBestScoringElementNoOptimal()
    {
        assertEquals(Optional.of(6), Utils.getBestScoringElement(Stream.of(1, 2, 6, 4, 5), s -> s, s -> s == 0, true));
        assertEquals(Optional.of(6), Utils.getBestScoringElement(Stream.of(6, 2, 1, 4, 5), s -> s, s -> s == 0, true));

        assertEquals(Optional.of(0), Utils.getBestScoringElement(Stream.of(5, 0, 3, 2, 1), s -> s, s -> s == 0, false));
        assertEquals(Optional.of(0), Utils.getBestScoringElement(Stream.of(0, 5, 3, 2, 1), s -> s, s -> s == 0, false));
    }
}
