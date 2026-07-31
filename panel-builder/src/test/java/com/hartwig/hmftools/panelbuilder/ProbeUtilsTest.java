package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.ProbeUtils.maxProbeEndOverlapping;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.maxProbeEndWithoutGap;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.minProbeStartOverlapping;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.minProbeStartWithoutGap;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionEndingAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionStartingAt;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class ProbeUtilsTest
{
    private static final int PROBE_START = 100;
    private static final int PROBE_CENTRE = 159;
    private static final int PROBE_END = 219;

    @Test
    public void testProbeRegionStartingAt()
    {
        assertEquals(new BaseRegion(PROBE_START, PROBE_END), probeRegionStartingAt(PROBE_START));
    }

    @Test
    public void testProbeRegionCenteredAt()
    {
        assertEquals(new BaseRegion(PROBE_START, PROBE_END), probeRegionCenteredAt(PROBE_CENTRE));
    }

    @Test
    public void testProbeRegionEndingAt()
    {
        assertEquals(new BaseRegion(PROBE_START, PROBE_END), probeRegionEndingAt(PROBE_END));
    }

    @Test
    public void testMinProbeStartContaining()
    {
        assertEquals(PROBE_START, ProbeUtils.minProbeStartContaining(PROBE_END));
    }

    @Test
    public void testMaxProbeEndContaining()
    {
        assertEquals(PROBE_END, ProbeUtils.maxProbeEndContaining(PROBE_START));
    }

    @Test
    public void testMinProbeStartOverlapping()
    {
        assertEquals(PROBE_START, minProbeStartOverlapping(new BaseRegion(PROBE_END, PROBE_END + 10)));
    }

    @Test
    public void testMaxProbeEndOverlapping()
    {
        assertEquals(PROBE_END, maxProbeEndOverlapping(new BaseRegion(PROBE_START - 10, PROBE_START)));
    }

    @Test
    public void testMinProbeStartWithoutGap()
    {
        assertEquals(PROBE_START, minProbeStartWithoutGap(new BaseRegion(PROBE_END + 1, PROBE_END + 10)));
    }

    @Test
    public void testMaxProbeEndWithoutGap()
    {
        assertEquals(PROBE_END, maxProbeEndWithoutGap(new BaseRegion(PROBE_START - 10, PROBE_START - 1)));
    }

    @Test
    public void testProbeTargetedRegionsSingleRegion()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.singleRegion(new ChrBaseRegion("1", 100, 199));
        TargetedRange targetedRange = new TargetedRange(10, 30);
        List<ChrBaseRegion> expected = List.of(new ChrBaseRegion("1", 110, 129));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsStartAndEndTargetedStart()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.variant(
                new ChrBaseRegion("1", 100, 199), Orientation.FORWARD,
                "",
                new ChrBaseRegion("2", 300, 399), Orientation.FORWARD);
        TargetedRange targetedRange = new TargetedRange(10, 30);
        List<ChrBaseRegion> expected = List.of(new ChrBaseRegion("1", 110, 129));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsStartAndEndTargetedEnd()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.variant(
                new ChrBaseRegion("1", 100, 199), Orientation.FORWARD,
                "",
                new ChrBaseRegion("2", 300, 399), Orientation.FORWARD);
        TargetedRange targetedRange = new TargetedRange(110, 130);
        List<ChrBaseRegion> expected = List.of(new ChrBaseRegion("2", 310, 329));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsStartAndEndTargetedBoth()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.variant(
                new ChrBaseRegion("1", 100, 199), Orientation.FORWARD,
                "",
                new ChrBaseRegion("2", 300, 399), Orientation.FORWARD);
        TargetedRange targetedRange = new TargetedRange(80, 130);
        List<ChrBaseRegion> expected = List.of(
                new ChrBaseRegion("1", 180, 199),
                new ChrBaseRegion("2", 300, 329));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsStartAndInsertTargetedStart()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.forwardSgl(
                new ChrBaseRegion("1", 100, 199), "AAAACCCCGGGGTTTTAAAA");
        TargetedRange targetedRange = new TargetedRange(10, 30);
        List<ChrBaseRegion> expected = List.of(new ChrBaseRegion("1", 110, 129));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsStartAndInsertTargetedInsert()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.forwardSgl(
                new ChrBaseRegion("1", 100, 199), "AAAACCCCGGGGTTTTAAAA");
        TargetedRange targetedRange = new TargetedRange(201, 210);
        List<ChrBaseRegion> expected = emptyList();
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsInsertAndEndTargetedInsert()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.reverseSgl(
                "AAAACCCCGGGGTTTTAAAA", new ChrBaseRegion("2", 300, 399));
        TargetedRange targetedRange = new TargetedRange(0, 16);
        List<ChrBaseRegion> expected = emptyList();
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsInsertAndEndTargetedEnd()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.reverseSgl(
                "AAAACCCCGGGGTTTTAAAA", new ChrBaseRegion("2", 300, 399));
        TargetedRange targetedRange = new TargetedRange(24, 54);
        List<ChrBaseRegion> expected = List.of(new ChrBaseRegion("2", 304, 333));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsStartInsertEndTargetedAll()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.variant(
                new ChrBaseRegion("1", 100, 199), Orientation.FORWARD,
                "AAAACCCCGGGGTTTTAAAA",
                new ChrBaseRegion("2", 300, 399), Orientation.FORWARD);
        TargetedRange targetedRange = new TargetedRange(20, 210);
        List<ChrBaseRegion> expected = List.of(
                new ChrBaseRegion("1", 120, 199),
                new ChrBaseRegion("2", 300, 389));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsStartInsertEndTargetedEnd()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.variant(
                new ChrBaseRegion("1", 100, 199), Orientation.FORWARD,
                "AAAACCCCGGGGTTTTAAAA",
                new ChrBaseRegion("2", 300, 399), Orientation.FORWARD);
        TargetedRange targetedRange = new TargetedRange(120, 150);
        List<ChrBaseRegion> expected = List.of(new ChrBaseRegion("2", 300, 329));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsThreeRegions()
    {
        SequenceDefinition sequenceDefinition = new SequenceDefinition(List.of(
                new RefSegment(new ChrBaseRegion("1", 100, 119), Orientation.FORWARD),
                new RefSegment(new ChrBaseRegion("2", 200, 219), Orientation.FORWARD),
                new RefSegment(new ChrBaseRegion("3", 300, 319), Orientation.FORWARD)));
        TargetedRange targetedRange = new TargetedRange(10, 50);
        List<ChrBaseRegion> expected = List.of(
                new ChrBaseRegion("1", 110, 119),
                new ChrBaseRegion("2", 200, 219),
                new ChrBaseRegion("3", 300, 309));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsReverseOrientation()
    {
        SequenceDefinition sequenceDefinition = SequenceDefinition.variant(
                new ChrBaseRegion("1", 100, 199), Orientation.REVERSE,
                "AAAACCCCGGGGTTTTAAAA",
                new ChrBaseRegion("2", 300, 399), Orientation.REVERSE);
        TargetedRange targetedRange = new TargetedRange(20, 210);
        List<ChrBaseRegion> expected = List.of(
                new ChrBaseRegion("1", 100, 179),
                new ChrBaseRegion("2", 310, 399));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsCrazyOverlapping()
    {
        // Genome-overlapping regions, mixed orientation, insert in the middle. The mapping is purely positional in probe-sequence space, so
        // genome overlap between segments is irrelevant and can even produce overlapping output regions.
        SequenceDefinition sequenceDefinition = new SequenceDefinition(List.of(
                new RefSegment(new ChrBaseRegion("1", 100, 109), Orientation.FORWARD),   // seq [0, 10)
                new InsertSeqSegment("ACGT"),                                            // seq [10, 14)
                new RefSegment(new ChrBaseRegion("1", 105, 114), Orientation.REVERSE),   // seq [14, 24), overlaps region 1
                new RefSegment(new ChrBaseRegion("2", 50, 54), Orientation.FORWARD)));   // seq [24, 29)
        TargetedRange targetedRange = new TargetedRange(5, 27);
        List<ChrBaseRegion> expected = List.of(
                // TODO: should collapse subsumed regions into 1?
                new ChrBaseRegion("1", 105, 109),   // region 1, offsets [5, 10)
                new ChrBaseRegion("1", 105, 114),   // region 2 reversed, offsets [0, 10)
                new ChrBaseRegion("2", 50, 52));    // region 3, offsets [0, 3)
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }
}
