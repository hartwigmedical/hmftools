package com.hartwig.hmftools.geneutils.paneldesign;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

// TODO

public class ProbeEvaluatorTest
{
    private static final ProbeEvaluator.Criteria CRITERIA = new ProbeEvaluator.Criteria(0.8, 0.45, 0.1);

    // TODO
    private static final ProbeQualityProfile PROBE_QUALITY_PROFILE = new ProbeQualityProfile(
            Map.of(),
            10, 5);

    private static final TargetRegion TARGET_REGION = new TargetRegion(
            new ChrBaseRegion("1", 100, 200),
            new TargetMetadata(TargetMetadata.Type.CUSTOM, "test"));

    private static final double EPSILON = 1e-9;

    private final ProbeEvaluator mProbeEvaluator;

    public ProbeEvaluatorTest()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(
                // 10b of good GC, 10b of good GC, 10b of low GC
                "1", "ACGTACGTACACGTACGTACAAAAAAAAAA"
        );
        // Compute chromosome lengths based on base sequences.
        refGenome.ChromosomeLengths.putAll(
                refGenome.RefGenomeMap.entrySet().stream().collect(
                        Collectors.toMap(Map.Entry::getKey, e -> e.getValue().length())));

        mProbeEvaluator = new ProbeEvaluator(refGenome, PROBE_QUALITY_PROFILE, null);
    }

    @Test
    public void testCriteriaGcContentMin()
    {
        assertEquals(0.35, CRITERIA.gcContentMin(), EPSILON);
    }

    @Test
    public void testCriteriaGcContentMax()
    {
        assertEquals(0.55, CRITERIA.gcContentMax(), EPSILON);
    }

    @Test
    public void testEvaluateProbesAcceptable()
    {
        Probe probe = new Probe(TARGET_REGION, new ChrBaseRegion("1", 1, 10));
        Probe evalProbe = mProbeEvaluator.evaluateProbe(probe, CRITERIA);
        assertNotSame(probe, evalProbe);
        assertEquals(probe.target(), evalProbe.target());
        assertEquals(CRITERIA, evalProbe.evalCriteria());
        assertNotNull(evalProbe.qualityScore());
        assertEquals(1, evalProbe.qualityScore(), EPSILON);
        assertEquals("ACGTACGTAC", evalProbe.sequence());
        assertNotNull(evalProbe.gcContent());
        assertEquals(0.5, evalProbe.gcContent(), EPSILON);
        assertNull(evalProbe.rejectionReason());
        assertTrue(evalProbe.accepted());
        assertFalse(evalProbe.rejected());
    }

    @Test
    public void testEvaluateProbesRejectQuality()
    {
        Probe probe = new Probe(TARGET_REGION, new ChrBaseRegion("1", 11, 20));
        Probe evalProbe = mProbeEvaluator.evaluateProbe(probe, CRITERIA);
        assertNotSame(probe, evalProbe);
        assertEquals(probe.target(), evalProbe.target());
        assertEquals(CRITERIA, evalProbe.evalCriteria());
        assertNotNull(evalProbe.qualityScore());
        assertEquals(0, evalProbe.qualityScore(), EPSILON);
        assertNotNull(evalProbe.rejectionReason());
        assertTrue(evalProbe.rejected());
        assertFalse(evalProbe.accepted());
        assertTrue(evalProbe.rejectionReason().toLowerCase().contains("quality"));
    }

    @Test
    public void testEvaluateProbesRejectGc()
    {
        Probe probe = new Probe(TARGET_REGION, new ChrBaseRegion("1", 21, 30));
        Probe evalProbe = mProbeEvaluator.evaluateProbe(probe, CRITERIA);
        assertNotSame(probe, evalProbe);
        assertEquals(probe.target(), evalProbe.target());
        assertEquals(CRITERIA, evalProbe.evalCriteria());
        assertEquals("AAAAAAAAAA", evalProbe.sequence());
        assertNotNull(evalProbe.gcContent());
        assertEquals(0, evalProbe.gcContent(), EPSILON);
        assertNotNull(evalProbe.rejectionReason());
        assertTrue(evalProbe.rejected());
        assertFalse(evalProbe.accepted());
        assertTrue(evalProbe.rejectionReason().toLowerCase().contains("gc"));
    }
}
