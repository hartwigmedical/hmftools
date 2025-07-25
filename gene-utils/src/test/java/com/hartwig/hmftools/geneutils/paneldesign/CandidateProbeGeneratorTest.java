package com.hartwig.hmftools.geneutils.paneldesign;

import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.probeRegionCenteredAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCentre;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class CandidateProbeGeneratorTest
{
    private static final TargetMetadata METADATA = new TargetMetadata(TargetMetadata.Type.CUSTOM, "extra");
    private static final ProbeContext CONTEXT = new ProbeContext(METADATA);

    private final CandidateProbeGenerator mGenerator = new CandidateProbeGenerator(Map.of(
            "1", 1_000_000,
            "2", PROBE_LENGTH - 1
    ));

    @Test
    public void testCoverOneSubregion()
    {
        ChrBaseRegion region = new ChrBaseRegion("1", 1000, 2000);

        List<Probe> actual = mGenerator.coverOneSubregion(region, CONTEXT).toList();

        Set<ChrBaseRegion> expectedRegions = IntStream.rangeClosed(region.start(), region.end() - PROBE_LENGTH + 1)
                .mapToObj(start -> new ChrBaseRegion(region.chromosome(), start, start + PROBE_LENGTH - 1))
                .collect(Collectors.toSet());
        List<ChrBaseRegion> actualRegions = actual.stream().map(Probe::region).toList();

        // Check the full set of probes is as expected.
        assertEquals(expectedRegions, new HashSet<>(actualRegions));

        // Check the first few are in the right order (no need to manually test all the rest).
        List<BaseRegion> actualFirstRegions = actualRegions.stream().map(ChrBaseRegion::baseRegion).limit(5).toList();
        int centre = regionCentre(region.baseRegion());
        List<BaseRegion> expectedFirstRegions = List.of(
                probeRegionCenteredAt(centre),
                probeRegionCenteredAt(centre + 1),
                probeRegionCenteredAt(centre - 1),
                probeRegionCenteredAt(centre + 2),
                probeRegionCenteredAt(centre - 2)
        );
        assertEquals(expectedFirstRegions, actualFirstRegions);

        assertProbeAttributes(actual);
    }

    @Test
    public void testCoverOneSubregionChromosomeBounds()
    {
        List<Probe> actual = mGenerator.coverOneSubregion(new ChrBaseRegion("1", 1, 1000), CONTEXT).toList();
        int minStart = actual.stream().mapToInt(probe -> probe.region().start()).min().orElseThrow();
        assertEquals(1, minStart);
    }

    @Test
    public void testCoverPosition()
    {
        String chromosome = "1";
        int position = 1000;
        Set<ChrBaseRegion> expectedRegions = IntStream.rangeClosed(position - PROBE_LENGTH + 1, position)
                .mapToObj(start -> new ChrBaseRegion(chromosome, start, start + PROBE_LENGTH - 1))
                .collect(Collectors.toSet());

        List<Probe> actual = mGenerator.coverPosition(new BasePosition(chromosome, position), CONTEXT).toList();
        List<ChrBaseRegion> actualRegions = actual.stream().map(Probe::region).toList();

        // Check the full set of probes is as expected.
        assertEquals(expectedRegions, new HashSet<>(actualRegions));

        // Check the first few are in the right order (no need to manually test all the rest).
        List<BaseRegion> actualFirstRegions = actualRegions.stream().map(ChrBaseRegion::baseRegion).limit(5).toList();
        List<BaseRegion> expectedFirstRegions = List.of(
                probeRegionCenteredAt(position),
                probeRegionCenteredAt(position + 1),
                probeRegionCenteredAt(position - 1),
                probeRegionCenteredAt(position + 2),
                probeRegionCenteredAt(position - 2)
        );
        assertEquals(expectedFirstRegions, actualFirstRegions);

        assertProbeAttributes(actual);
    }

    @Test
    public void testCoverPositionChromosomeBounds()
    {
        List<Probe> actual = mGenerator.coverPosition(new BasePosition("2", 1), CONTEXT).toList();
        assertEquals(emptyList(), actual);
    }

    @Test
    public void testAllOverlapping()
    {
        String chromosome = "1";
        ChrBaseRegion region = new ChrBaseRegion(chromosome, 1000, 2000);
        List<ChrBaseRegion> expectedRegions = IntStream.rangeClosed(region.start() - PROBE_LENGTH + 1, region.end())
                .mapToObj(start -> new ChrBaseRegion(chromosome, start, start + PROBE_LENGTH - 1))
                .toList();

        List<Probe> actual = mGenerator.allOverlapping(region, CONTEXT).toList();
        List<ChrBaseRegion> actualRegions = actual.stream().map(Probe::region).toList();

        assertEquals(expectedRegions, actualRegions);

        assertProbeAttributes(actual);
    }

    @Test
    public void testAllOverlappingChromosomeBounds()
    {
        List<Probe> actual = mGenerator.allOverlapping(new ChrBaseRegion("2", 1, 10), CONTEXT).toList();
        assertEquals(emptyList(), actual);
    }

    private static void assertProbeAttributes(final List<Probe> probes)
    {
        probes.forEach(probe ->
        {
            assertEquals(METADATA, probe.metadata());
            assertFalse(probe.evaluated());
        });
    }
}
