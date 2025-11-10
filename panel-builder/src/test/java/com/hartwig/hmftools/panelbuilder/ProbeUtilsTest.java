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
        SequenceDefinition sequenceDefinition = new SequenceDefinition(
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
        SequenceDefinition sequenceDefinition = new SequenceDefinition(
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
        SequenceDefinition sequenceDefinition = new SequenceDefinition(
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
        SequenceDefinition sequenceDefinition = new SequenceDefinition(
                new ChrBaseRegion("1", 100, 199), Orientation.FORWARD,
                "AAAACCCCGGGGTTTTAAAA",
                null, null);
        TargetedRange targetedRange = new TargetedRange(10, 30);
        List<ChrBaseRegion> expected = List.of(new ChrBaseRegion("1", 110, 129));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsStartAndInsertTargetedInsert()
    {
        SequenceDefinition sequenceDefinition = new SequenceDefinition(
                new ChrBaseRegion("1", 100, 199), Orientation.FORWARD,
                "AAAACCCCGGGGTTTTAAAA",
                null, null);
        TargetedRange targetedRange = new TargetedRange(201, 210);
        List<ChrBaseRegion> expected = emptyList();
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsInsertAndEndTargetedInsert()
    {
        SequenceDefinition sequenceDefinition = new SequenceDefinition(
                null, null,
                "AAAACCCCGGGGTTTTAAAA",
                new ChrBaseRegion("2", 300, 399), Orientation.FORWARD);
        TargetedRange targetedRange = new TargetedRange(0, 16);
        List<ChrBaseRegion> expected = emptyList();
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsInsertAndEndTargetedEnd()
    {
        SequenceDefinition sequenceDefinition = new SequenceDefinition(
                null, null,
                "AAAACCCCGGGGTTTTAAAA",
                new ChrBaseRegion("2", 300, 399), Orientation.FORWARD);
        TargetedRange targetedRange = new TargetedRange(24, 54);
        List<ChrBaseRegion> expected = List.of(new ChrBaseRegion("2", 304, 333));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsStartInsertEndTargetedAll()
    {
        SequenceDefinition sequenceDefinition = new SequenceDefinition(
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
        SequenceDefinition sequenceDefinition = new SequenceDefinition(
                new ChrBaseRegion("1", 100, 199), Orientation.FORWARD,
                "AAAACCCCGGGGTTTTAAAA",
                new ChrBaseRegion("2", 300, 399), Orientation.FORWARD);
        TargetedRange targetedRange = new TargetedRange(120, 150);
        List<ChrBaseRegion> expected = List.of(new ChrBaseRegion("2", 300, 329));
        List<ChrBaseRegion> actual = ProbeUtils.probeTargetedRegions(sequenceDefinition, targetedRange);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeTargetedRegionsReverseOrientation()
    {
        SequenceDefinition sequenceDefinition = new SequenceDefinition(
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
}
