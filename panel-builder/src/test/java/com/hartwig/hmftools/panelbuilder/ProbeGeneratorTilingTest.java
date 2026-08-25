package com.hartwig.hmftools.panelbuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.IntPredicate;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

// Characterisation tests for the DNA whole-region tiling path (ProbeGenerator.coverRegion -> coverUncoveredRegion ->
// coverAcceptableSubregion). These lock in current behaviour so the planned merge onto a shared mapping-driven generator can be validated.
// PROBE_LENGTH is 120, REGION_UNCOVERED_MAX is 20, PROBE_OVERLAP_EXTENSION_BALANCE is 0.5, PROBE_SHIFT_MAX is 5.
public class ProbeGeneratorTilingTest
{
    private static final int CHR_LENGTH = 10_000;
    private static final TargetMetadata METADATA = new TargetMetadata(TargetMetadata.Type.CUSTOM_REGION, "dna");
    // Quality score alone decides acceptance; GC is always within tolerance.
    private static final ProbeEvaluator.Criteria CRITERIA = new ProbeEvaluator.Criteria(0.5, 0.5, 1.0);
    private static final ProbeSelector.Strategy SELECT = new ProbeSelector.Strategy.FirstAcceptable();

    @Test
    public void testExactMultipleAllAcceptable()
    {
        // Region length an exact multiple of the probe length: probes tile flush with no overlap or extension.
        ProbeGenerationResult result = generate(new ChrBaseRegion("1", 1001, 1360), start -> true);

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
        // Region a little longer than a probe: one probe centred on the region, small uncovered edges permitted (not rejected).
        ProbeGenerationResult result = generate(new ChrBaseRegion("1", 1001, 1130), start -> true);

        assertEquals(List.of(new ChrBaseRegion("1", 1006, 1125)), singleRegions(result));
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testMultipleProbesCentredWithEdgeSlack()
    {
        // Region just over twice the probe length: two probes centred with a few uncovered bases at each edge (NOT flushed to the edges;
        // this is the key DNA behaviour vs RNA edge-pinning).
        ProbeGenerationResult result = generate(new ChrBaseRegion("1", 1001, 1250), start -> true);

        assertEquals(
                List.of(new ChrBaseRegion("1", 1006, 1125), new ChrBaseRegion("1", 1126, 1245)),
                singleRegions(result));
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testProbesExtendBeyondRegion()
    {
        // Region length leaves "extra" probe bases which are split between overlap and extension past the region edges (governed by
        // PROBE_OVERLAP_EXTENSION_BALANCE). Probes here extend slightly outside [1001, 1290].
        ProbeGenerationResult result = generate(new ChrBaseRegion("1", 1001, 1290), start -> true);

        assertEquals(
                List.of(
                        new ChrBaseRegion("1", 983, 1102),
                        new ChrBaseRegion("1", 1086, 1205),
                        new ChrBaseRegion("1", 1188, 1307)),
                singleRegions(result));
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testRegionSmallerThanProbe()
    {
        // Region shorter than a probe: a single probe centred on it, extending beyond both edges (a DNA probe can capture flanking sequence,
        // unlike RNA which must stay within exons).
        ProbeGenerationResult result = generate(new ChrBaseRegion("1", 1001, 1050), start -> true);

        assertEquals(List.of(new ChrBaseRegion("1", 966, 1085)), singleRegions(result));
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    @Test
    public void testRejectedRegionSplit()
    {
        // A rejected band in the middle splits the region into two acceptable subregions, each tiled separately (centred with extension,
        // clamped to the subregion); the band is rejected. Reject probes overlapping genome [1191, 1211], i.e. starts in [1072, 1211].
        // Each subregion tiles two probes: the probe abutting the rejection is flush to it (1190, 1212), the outer probe extends past the
        // target edge (to 976 and 1426) - the DNA "probes may extend beyond the region" behaviour.
        IntPredicate acceptable = start -> start <= 1071 || start >= 1212;
        ProbeGenerationResult result = generate(new ChrBaseRegion("1", 1001, 1400), acceptable);

        assertEquals(
                List.of(
                        new ChrBaseRegion("1", 976, 1095),
                        new ChrBaseRegion("1", 1071, 1190),
                        new ChrBaseRegion("1", 1212, 1331),
                        new ChrBaseRegion("1", 1307, 1426)),
                singleRegions(result));
        assertEquals(
                List.of(new ChrBaseRegion("1", 1191, 1211)),
                result.rejectedFeatures().stream().map(RejectedFeature::region).toList());
    }

    @Test
    public void testWholeRegionRejected()
    {
        // No acceptable probe anywhere: no probes, whole region rejected.
        ProbeGenerationResult result = generate(new ChrBaseRegion("1", 1001, 1360), start -> false);

        assertTrue(result.probes().isEmpty());
        assertEquals(
                List.of(new ChrBaseRegion("1", 1001, 1360)),
                result.rejectedFeatures().stream().map(RejectedFeature::region).toList());
    }

    @Test
    public void testChromosomeStartClamp()
    {
        // Region at the chromosome start: probes cannot start before position 1.
        ProbeGenerationResult result = generate(new ChrBaseRegion("1", 1, 200), start -> true);

        List<ChrBaseRegion> probes = singleRegions(result);
        assertEquals(1, probes.get(0).start());
        assertTrue(result.rejectedFeatures().isEmpty());
    }

    private static ProbeGenerationResult generate(final ChrBaseRegion region, final IntPredicate acceptableStart)
    {
        ProbeGenerator generator = new ProbeGenerator(Map.of("1", CHR_LENGTH), fakeEvaluator(acceptableStart), null);
        ProbeGenerationSpec spec = new ProbeGenerationSpec.CoverRegion(region, METADATA, CRITERIA, SELECT);
        return generator.generateBatch(Stream.of(spec), new PanelData());
    }

    // Evaluator with controllable acceptance: a probe is accepted iff its genome start satisfies the predicate (via quality score). GC is
    // always 0.5 (within tolerance), so quality alone decides.
    private static ProbeEvaluator fakeEvaluator(final IntPredicate acceptableStart)
    {
        return new ProbeEvaluator(
                probe ->
                {
                    byte[] bases = "A".repeat(probe.definition().baseLength()).getBytes();
                    return new SequenceData(bases, true, bases.length / 2);
                },
                probes -> probes.map(probe ->
                        probe.withQualityScore(acceptableStart.test(probe.definition().singleRegion().start()) ? 1.0 : 0.0)));
    }

    private static List<ChrBaseRegion> singleRegions(final ProbeGenerationResult result)
    {
        return result.probes().stream()
                .map(probe -> probe.definition().singleRegion())
                .sorted(Comparator.comparingInt(ChrBaseRegion::start))
                .toList();
    }
}
