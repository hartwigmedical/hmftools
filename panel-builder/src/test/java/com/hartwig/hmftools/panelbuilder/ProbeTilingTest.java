package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.ProbeTiling.calculateContainedTiling;
import static com.hartwig.hmftools.panelbuilder.ProbeTiling.calculateOptimalProbeTiling;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import org.junit.Test;

// PROBE_LENGTH is 120 and REGION_UNCOVERED_MAX is 20 (see PanelBuilderConstants).
public class ProbeTilingTest
{
    @Test
    public void testPerfectSingleProbe()
    {
        // Region exactly one probe long: one probe starting at the region start.
        assertEquals(List.of(1), calculateOptimalProbeTiling(new BaseRegion(1, 120), new BaseRegion(1, 120)));
    }

    @Test
    public void testShorterThanProbeCentred()
    {
        // Region shorter than a probe: one probe, centred on the region within the bounds.
        assertEquals(List.of(1), calculateOptimalProbeTiling(new BaseRegion(1, 100), new BaseRegion(1, 300)));
    }

    @Test
    public void testTwoProbes()
    {
        // Region 200 long: two probes, overlap and extension balanced.
        assertEquals(List.of(1, 101), calculateOptimalProbeTiling(new BaseRegion(1, 200), new BaseRegion(1, 400)));
    }

    @Test
    public void testNoTilingPossibleWhenBoundsTooTight()
    {
        // Region far shorter than a probe and bounds leave no room to place a full probe.
        assertTrue(calculateOptimalProbeTiling(new BaseRegion(1, 50), new BaseRegion(1, 50)).isEmpty());
    }

    @Test
    public void testProbeBoundsMustContainRegion()
    {
        assertThrows(
                IllegalArgumentException.class,
                () -> calculateOptimalProbeTiling(new BaseRegion(1, 200), new BaseRegion(1, 100)));
    }

    @Test
    public void testContainedBothPinned()
    {
        // Both edges pinned flush: first probe starts at region start, last probe ends at region end.
        assertEquals(List.of(1, 81), calculateContainedTiling(new BaseRegion(1, 200), true, true));
        // Perfect tiling: three abutting probes flush to both ends.
        assertEquals(List.of(1, 121, 241), calculateContainedTiling(new BaseRegion(1, 360), true, true));
    }

    @Test
    public void testContainedPinStartOnly()
    {
        // Pinned left, small uncovered gap allowed at the unpinned right edge.
        assertEquals(List.of(1), calculateContainedTiling(new BaseRegion(1, 130), true, false));
    }

    @Test
    public void testContainedPinEndOnly()
    {
        // Pinned right (last probe ends at region end), gap allowed at the unpinned left edge.
        assertEquals(List.of(11), calculateContainedTiling(new BaseRegion(1, 130), false, true));
    }

    @Test
    public void testContainedNeitherPinned()
    {
        // Neither edge pinned (e.g. a sub-range between two rejected regions): centred, gap split between edges.
        assertEquals(List.of(6), calculateContainedTiling(new BaseRegion(1, 130), false, false));
    }

    @Test
    public void testContainedTooShortToTile()
    {
        // Shorter than a probe: no tiling (caller pads instead).
        assertTrue(calculateContainedTiling(new BaseRegion(1, 100), true, true).isEmpty());
    }
}
