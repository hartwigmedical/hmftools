package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCentre;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class CandidateProbeGeneratorTest
{
    private static final TargetMetadata METADATA = new TargetMetadata(TargetMetadata.Type.CUSTOM, "extra");

    private final CandidateProbeGenerator mGenerator;

    public CandidateProbeGeneratorTest()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.ChromosomeLengths = Map.of(
                "1", 10_000,
                "2", PROBE_LENGTH - 1
        );
        for(Map.Entry<String, Integer> entry : refGenome.ChromosomeLengths.entrySet())
        {
            refGenome.RefGenomeMap.put(entry.getKey(), MockRefGenome.generateRandomBases(entry.getValue()));
        }
        mGenerator = new CandidateProbeGenerator(refGenome.chromosomeLengths());
    }

    @Test
    public void testCoverOneSubregion()
    {
        ChrBaseRegion region = new ChrBaseRegion("1", 1000, 2000);

        List<Probe> actual = mGenerator.coverOneSubregion(region, METADATA).toList();

        Set<ChrBaseRegion> expectedRegions = IntStream.rangeClosed(region.start(), region.end() - PROBE_LENGTH + 1)
                .mapToObj(start -> new ChrBaseRegion(region.chromosome(), start, start + PROBE_LENGTH - 1))
                .collect(Collectors.toSet());
        List<ChrBaseRegion> actualRegions = actual.stream().map(probe -> probe.definition().singleRegion()).toList();

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
        List<Probe> actual = mGenerator.coverOneSubregion(new ChrBaseRegion("1", 1, 1000), METADATA).toList();
        int minStart = actual.stream().mapToInt(probe -> probe.definition().singleRegion().start()).min().orElseThrow();
        assertEquals(1, minStart);
    }

    @Test
    public void testAllOverlapping()
    {
        String chromosome = "1";
        ChrBaseRegion region = new ChrBaseRegion(chromosome, 1000, 2000);
        List<ChrBaseRegion> expectedRegions = IntStream.rangeClosed(region.start() - PROBE_LENGTH + 1, region.end())
                .mapToObj(start -> new ChrBaseRegion(chromosome, start, start + PROBE_LENGTH - 1))
                .toList();

        List<Probe> actual = mGenerator.allOverlapping(region, METADATA).toList();
        List<ChrBaseRegion> actualRegions = actual.stream().map(probe -> probe.definition().singleRegion()).toList();

        assertEquals(expectedRegions, actualRegions);

        assertProbeAttributes(actual);
    }

    @Test
    public void testAllOverlappingChromosomeBounds()
    {
        List<Probe> actual = mGenerator.allOverlapping(new ChrBaseRegion("2", 1, 10), METADATA).toList();
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
