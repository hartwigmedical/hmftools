package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.function.IntPredicate;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

// Exon-aware (RNA) tiling via ProbeGenerator.coverExonRange: probe-space tiling over an exon RegionMapping, pinned flush to exon boundaries,
// spliced probes across junctions.
// PROBE_LENGTH is 120, REGION_UNCOVERED_MAX is 20, RNA_EXON_SINGLE_PROBE_SLACK is 20, PROBE_SHIFT_MAX is 5.
// Exon genome coordinates are chosen so that probe-space position s maps to genome position 1001+s within the first exon.
public class ProbeGeneratorExonRangeTest
{
    private static final TargetMetadata METADATA = new TargetMetadata(TargetMetadata.Type.CUSTOM_REGION, "rna");
    // Quality score alone decides acceptance; GC is always within tolerance.
    private static final ProbeEvaluator.Criteria CRITERIA = new ProbeEvaluator.Criteria(0.5, 0.5, 1.0);
    private static final ProbeSelector.Strategy STRATEGY = new ProbeSelector.Strategy.FirstAcceptable();

    @Test
    public void testMultiProbePinnedTiling()
    {
        // Whole single exon, fully acceptable, length a multiple of the probe length: probes tiled flush to both exon boundaries.
        RegionMapping mapping = new RegionMapping(List.of(new ChrBaseRegion("1", 1001, 1360)));
        ProbeGenerationResult result = generate(mapping, 0, 360, start -> true);

        assertEquals(
                List.of(
                        new ChrBaseRegion("1", 1001, 1120),
                        new ChrBaseRegion("1", 1121, 1240),
                        new ChrBaseRegion("1", 1241, 1360)),
                singleRegions(result));
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testSingleCentredProbe()
    {
        // Exon between the probe length and probe length + slack: a single probe centred in the exon.
        RegionMapping mapping = new RegionMapping(List.of(new ChrBaseRegion("1", 1001, 1130)));
        ProbeGenerationResult result = generate(mapping, 0, 130, start -> true);

        // Centred: (130 - 120) / 2 = 5 offset, so probe-space [5, 125) -> genome [1006, 1125].
        assertEquals(List.of(new ChrBaseRegion("1", 1006, 1125)), singleRegions(result));
        // Edge slack at the unpinned edges is permitted, not rejected.
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testUnpinnedInteriorTiling()
    {
        // Target range in the middle of a large exon (neither edge is an exon boundary): tiling centred within the range with edge slack.
        RegionMapping mapping = new RegionMapping(List.of(new ChrBaseRegion("1", 1001, 1600)));
        ProbeGenerationResult result = generate(mapping, 100, 460, start -> true);

        assertEquals(
                List.of(
                        new ChrBaseRegion("1", 1101, 1220),
                        new ChrBaseRegion("1", 1221, 1340),
                        new ChrBaseRegion("1", 1341, 1460)),
                singleRegions(result));
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testUnpinnedInteriorRangeUncoveredBudget()
    {
        // Interior target range (both edges unpinned), length 145. REGION_UNCOVERED_MAX (20) is the total uncovered budget, so two probes are
        // needed. With both edges interior (not exon boundaries), the probes are centred and may extend past the target into the surrounding
        // exon sequence (the same behaviour as DNA tiling). GenesRna only targets whole exons, so this interior-target case is synthetic.
        RegionMapping mapping = new RegionMapping(List.of(new ChrBaseRegion("1", 1001, 1600)));
        ProbeGenerationResult result = generate(mapping, 200, 345, start -> true);

        assertEquals(
                List.of(new ChrBaseRegion("1", 1177, 1296), new ChrBaseRegion("1", 1249, 1368)),
                singleRegions(result));
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testShortExonBothNeighbours()
    {
        // Short exon with an adjacent exon on both sides: one centred spliced probe pulls padding evenly from both neighbours.
        RegionMapping mapping = new RegionMapping(List.of(
                new ChrBaseRegion("1", 1001, 1200),
                new ChrBaseRegion("1", 2001, 2060),
                new ChrBaseRegion("1", 3001, 3200)));
        // Target the short middle exon, probe-space [200, 260), length 60. 60 padding needed, split 30/30 across the two junctions.
        ProbeGenerationResult result = generate(mapping, 200, 260, start -> true);

        assertEquals(1, result.probes().size());
        assertEquals(
                List.of(
                        new ChrBaseRegion("1", 1171, 1200),
                        new ChrBaseRegion("1", 2001, 2060),
                        new ChrBaseRegion("1", 3001, 3030)),
                result.probes().get(0).definition().regions());
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testShortExonLeftNeighbourOnly()
    {
        // Short exon at the transcript end (only a left neighbour): padding is asymmetric, all pulled from the preceding exon.
        RegionMapping mapping = new RegionMapping(List.of(
                new ChrBaseRegion("1", 1001, 1200),
                new ChrBaseRegion("1", 2001, 2060)));
        // Target the short last exon, probe-space [200, 260), length 60. All 60 padding comes from the left exon.
        ProbeGenerationResult result = generate(mapping, 200, 260, start -> true);

        assertEquals(1, result.probes().size());
        assertEquals(
                List.of(new ChrBaseRegion("1", 1141, 1200), new ChrBaseRegion("1", 2001, 2060)),
                result.probes().get(0).definition().regions());
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testShortExonRightNeighbourOnly()
    {
        // Short exon at the transcript start (only a right neighbour): padding is asymmetric, all pulled from the following exon.
        RegionMapping mapping = new RegionMapping(List.of(
                new ChrBaseRegion("1", 1001, 1060),
                new ChrBaseRegion("1", 2001, 2200)));
        // Target the short first exon, probe-space [0, 60), length 60. All 60 padding comes from the right exon.
        ProbeGenerationResult result = generate(mapping, 0, 60, start -> true);

        assertEquals(1, result.probes().size());
        assertEquals(
                List.of(new ChrBaseRegion("1", 1001, 1060), new ChrBaseRegion("1", 2001, 2060)),
                result.probes().get(0).definition().regions());
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testShortExonShortNeighbourAtTerminus()
    {
        // Short target exon whose only right neighbour is itself short and is the last exon (terminus beyond it). The right side cannot
        // supply the full padding, so the probe shifts left and takes the deficit from the larger left exon (asymmetric).
        RegionMapping mapping = new RegionMapping(List.of(
                new ChrBaseRegion("1", 1001, 1200),
                new ChrBaseRegion("1", 2001, 2040),
                new ChrBaseRegion("1", 3001, 3020)));
        // Target the middle exon, probe-space [200, 240), length 40. Needs 80 padding; right exon only has 20, so 60 comes from the left.
        ProbeGenerationResult result = generate(mapping, 200, 240, start -> true);

        assertEquals(1, result.probes().size());
        assertEquals(
                List.of(
                        new ChrBaseRegion("1", 1141, 1200),
                        new ChrBaseRegion("1", 2001, 2040),
                        new ChrBaseRegion("1", 3001, 3020)),
                result.probes().get(0).definition().regions());
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testShortExonShortNeighbourNotTerminus()
    {
        // Short target exon whose right neighbour is short but is followed by a further exon. Padding stays symmetric because the right side
        // has enough sequence overall, spanning the short neighbour and into the second-adjacent exon.
        RegionMapping mapping = new RegionMapping(List.of(
                new ChrBaseRegion("1", 1001, 1200),
                new ChrBaseRegion("1", 2001, 2040),
                new ChrBaseRegion("1", 3001, 3010),
                new ChrBaseRegion("1", 4001, 4200)));
        // Target the second exon, probe-space [200, 240), length 40. Symmetric 40/40 padding; the right 40 spans the 10b exon then 30 more.
        ProbeGenerationResult result = generate(mapping, 200, 240, start -> true);

        assertEquals(1, result.probes().size());
        assertEquals(
                List.of(
                        new ChrBaseRegion("1", 1161, 1200),
                        new ChrBaseRegion("1", 2001, 2040),
                        new ChrBaseRegion("1", 3001, 3010),
                        new ChrBaseRegion("1", 4001, 4030)),
                result.probes().get(0).definition().regions());
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testShortExonBothNeighboursTooShort()
    {
        // Short target exon flanked by short terminal exons: the whole transcript is shorter than a probe, so no padding probe can be placed.
        RegionMapping mapping = new RegionMapping(List.of(
                new ChrBaseRegion("1", 1001, 1030),
                new ChrBaseRegion("1", 2001, 2040),
                new ChrBaseRegion("1", 3001, 3030)));
        // Target the middle exon, probe-space [30, 70), length 40. Total mapped length is only 100 (< probe length), so no probe fits.
        ProbeGenerationResult result = generate(mapping, 30, 70, start -> true);

        assertTrue(result.probes().isEmpty());
        assertEquals(
                List.of(new ChrBaseRegion("1", 2001, 2040)),
                result.rejectedFeatures().stream().map(RejectedFeature::region).toList());
    }

    @Test
    public void testStrandedShortPieceNoNeighbour()
    {
        // A rejected band strands a short piece against the transcript terminus. With no adjacent exon to pad from, the piece cannot form an
        // acceptable (spliced) sub-range at all, so it is simply rejected - no probe is attempted off the end of the mapping.
        RegionMapping mapping = new RegionMapping(List.of(new ChrBaseRegion("1", 1001, 1300)));
        // Reject everything from probe-space start 111 onward, leaving acceptable windows only in [0, 110].
        ProbeGenerationResult result = generate(mapping, 0, 300, start -> start <= 110);

        assertEquals(
                List.of(new ChrBaseRegion("1", 1001, 1120), new ChrBaseRegion("1", 1111, 1230)),
                singleRegions(result));
        // The stranded tail [1231, 1300] is rejected rather than covered by a probe running off the end.
        assertEquals(
                List.of(new ChrBaseRegion("1", 1231, 1300)),
                result.rejectedFeatures().stream().map(RejectedFeature::region).toList());
    }

    @Test
    public void testExonJustOverProbeMultipleFlushEdges()
    {
        // Exon a bit longer than twice the probe length, both edges exon boundaries: three probes, outermost flush to the exon edges.
        RegionMapping mapping = new RegionMapping(List.of(new ChrBaseRegion("1", 1001, 1250)));
        ProbeGenerationResult result = generate(mapping, 0, 250, start -> true);

        List<ChrBaseRegion> probes = singleRegions(result);
        assertEquals(
                List.of(
                        new ChrBaseRegion("1", 1001, 1120),
                        new ChrBaseRegion("1", 1066, 1185),
                        new ChrBaseRegion("1", 1131, 1250)),
                probes);
        // First probe starts at the exon start, last probe ends at the exon end.
        assertEquals(1001, probes.get(0).start());
        assertEquals(1250, probes.get(probes.size() - 1).end());
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testInvalidRangeOutsideMapping()
    {
        // Range extending beyond the mapped length.
        RegionMapping singleExon = new RegionMapping(List.of(new ChrBaseRegion("1", 1001, 1300)));
        assertThrows(IllegalArgumentException.class, () -> generate(singleExon, 0, 301, start -> true));

        // Range spanning a junction (not within a single exon).
        RegionMapping twoExons = new RegionMapping(List.of(
                new ChrBaseRegion("1", 1001, 1200),
                new ChrBaseRegion("1", 2001, 2200)));
        assertThrows(IllegalArgumentException.class, () -> generate(twoExons, 150, 250, start -> true));
    }

    @Test
    public void testRejectedRegionSplit()
    {
        // A rejected band in the middle of an exon splits it into two acceptable sub-ranges, each tiled separately; the band is rejected.
        RegionMapping mapping = new RegionMapping(List.of(new ChrBaseRegion("1", 1001, 1400)));
        // Reject windows overlapping probe-space [190, 210], i.e. starts in [71, 210].
        IntPredicate acceptable = start -> start <= 70 || start >= 211;
        ProbeGenerationResult result = generate(mapping, 0, 400, acceptable);

        assertEquals(
                List.of(
                        new ChrBaseRegion("1", 1001, 1120),
                        new ChrBaseRegion("1", 1071, 1190),
                        new ChrBaseRegion("1", 1212, 1331),
                        new ChrBaseRegion("1", 1281, 1400)),
                singleRegions(result));
        // The uncovered band between the two acceptable sub-ranges is rejected.
        assertEquals(
                List.of(new ChrBaseRegion("1", 1191, 1211)),
                result.rejectedFeatures().stream().map(RejectedFeature::region).toList());
    }

    @Test
    public void testWholeRangeRejected()
    {
        // No acceptable probe anywhere: no probes, whole target rejected.
        RegionMapping mapping = new RegionMapping(List.of(new ChrBaseRegion("1", 1001, 1360)));
        ProbeGenerationResult result = generate(mapping, 0, 360, start -> false);

        assertTrue(result.probes().isEmpty());
        assertEquals(
                List.of(new ChrBaseRegion("1", 1001, 1360)),
                result.rejectedFeatures().stream().map(RejectedFeature::region).toList());
    }

    // Exercises the batch (production) interface with a single spec; the single-spec convenience lives here in test code, not production.
    private static ProbeGenerationResult generate(final RegionMapping mapping, int rangeStart, int rangeEnd,
            final IntPredicate acceptableStart)
    {
        ProbeGenerator generator = new ProbeGenerator(Map.of(), fakeEvaluator(mapping, acceptableStart), null);
        ProbeGenerationSpec spec = new ProbeGenerationSpec.CoverExonRange(mapping, rangeStart, rangeEnd, METADATA, CRITERIA, STRATEGY);
        ProbeGenerationResult result = generator.generateBatch(Stream.of(spec), new NoCoveragePanel());
        assertProbeInvariants(mapping, rangeStart, rangeEnd, result);
        return result;
    }

    // Minimal PanelBuffer: RNA exon-range generation does no coverage subtraction, and the result is taken from generateBatch's return value.
    private static class NoCoveragePanel implements PanelBuffer
    {
        @Override
        public Stream<ChrBaseRegion> coveredRegions()
        {
            return Stream.empty();
        }

        @Override
        public void addResult(final ProbeGenerationResult result)
        {
        }
    }

    // Generic invariants that must hold for any result, checked automatically for every scenario.
    private static void assertProbeInvariants(final RegionMapping mapping, int rangeStart, int rangeEnd, final ProbeGenerationResult result)
    {
        List<Integer> starts = new ArrayList<>();
        for(Probe probe : result.probes())
        {
            // Each probe covers exactly a probe-length window fully contained in the mapping, its regions being the exon-boundary split of
            // that window (so it never covers intronic / off-transcript bases).
            int start = probeSpaceStart(mapping, probe);
            starts.add(start);
            assertEquals(PROBE_LENGTH, probe.definition().baseLength());
            assertTrue(start >= 0 && start + PROBE_LENGTH <= mapping.length());
            assertEquals(mapping.toGenomeRegions(start, start + PROBE_LENGTH), probe.definition().regions());
            // The probe only targets bases within the requested range.
            assertTrue(start + probe.targetedRange().startOffset() >= rangeStart);
            assertTrue(start + probe.targetedRange().endOffset() <= rangeEnd);
            // For a spliced probe, its regions are strictly ascending and non-overlapping (verified independently of toGenomeRegions).
            List<ChrBaseRegion> regions = probe.definition().regions();
            for(int j = 1; j < regions.size(); ++j)
            {
                assertTrue(regions.get(j - 1).compareTo(regions.get(j)) < 0 && !regions.get(j - 1).overlaps(regions.get(j)));
            }
        }

        // Probes are ordered by probe-space start with no duplicates, so consecutive probes overlap by less than a probe length.
        assertEquals(starts.stream().sorted().distinct().toList(), starts);

        // Rejected regions lie within the mapping (never intronic / off-transcript); collect them in probe-space for the gap check.
        List<BaseRegion> rejectedRanges = new ArrayList<>();
        for(RejectedFeature rejected : result.rejectedFeatures())
        {
            ChrBaseRegion region = rejected.region();
            OptionalInt spaceStart = mapping.toProbeSpacePosition(region.chromosome(), region.start());
            OptionalInt spaceEnd = mapping.toProbeSpacePosition(region.chromosome(), region.end());
            assertTrue(spaceStart.isPresent() && spaceEnd.isPresent());
            rejectedRanges.add(new BaseRegion(spaceStart.getAsInt(), spaceEnd.getAsInt()));
        }

        // Within a tiling there are never gaps between probes, so any gap between consecutive probes must be an entirely rejected region.
        for(int i = 1; i < starts.size(); ++i)
        {
            int gapStart = starts.get(i - 1) + PROBE_LENGTH;
            int gapEnd = starts.get(i) - 1;
            if(gapStart <= gapEnd)
            {
                assertTrue(rejectedRanges.stream().anyMatch(range -> range.start() <= gapStart && range.end() >= gapEnd));
            }
        }
    }

    private static int probeSpaceStart(final RegionMapping mapping, final Probe probe)
    {
        ChrBaseRegion first = probe.definition().regions().get(0);
        return mapping.toProbeSpacePosition(first.chromosome(), first.start()).orElseThrow();
    }

    // Evaluator with controllable acceptance: a probe is accepted iff its probe-space start satisfies the predicate (via quality score).
    private static ProbeEvaluator fakeEvaluator(final RegionMapping mapping, final IntPredicate acceptableStart)
    {
        return new ProbeEvaluator(
                probe ->
                {
                    // Fake sequence and GC - doesn't matter,
                    byte[] bases = "A".repeat(probe.definition().baseLength()).getBytes();
                    return new SequenceData(bases, true, bases.length / 2);
                },
                probes -> probes.map(probe ->
                {
                    // Accept or reject the probe based on its start position, as specified by the caller.
                    return probe.withQualityScore(acceptableStart.test(probeSpaceStart(mapping, probe)) ? 1.0 : 0.0);
                }));
    }

    private static List<ChrBaseRegion> singleRegions(final ProbeGenerationResult result)
    {
        return result.probes().stream()
                .map(probe -> probe.definition().singleRegion())
                .sorted(Comparator.comparingInt(ChrBaseRegion::start))
                .toList();
    }
}
