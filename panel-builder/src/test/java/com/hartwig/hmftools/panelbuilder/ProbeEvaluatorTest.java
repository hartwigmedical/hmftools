package com.hartwig.hmftools.panelbuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.mappability.ProbeQualityWindow;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class ProbeEvaluatorTest
{
    private static final ProbeEvaluator.Criteria CRITERIA = new ProbeEvaluator.Criteria(0.8, 0.45, 0.1);

    private static final TargetMetadata METADATA = new TargetMetadata(TargetMetadata.Type.CUSTOM, "test");

    private static final double EPSILON = 1e-6;

    private final ProbeEvaluator mProbeEvaluator;

    public ProbeEvaluatorTest()
    {
        // Genome setup:
        //   1-10: acceptable quality and CG
        //   11-20: acceptable quality, rejected GC
        //   21-30: rejected quality, accepted GC

        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(
                // 10b of acceptable GC, 10b of rejected GC, 10b of accepted GC
                "1", "ACGTACGTACAAAAAAAAAAACGTACGTAC"
        );
        // Compute chromosome lengths based on base sequences.
        refGenome.ChromosomeLengths.putAll(
                refGenome.RefGenomeMap.entrySet().stream().collect(
                        Collectors.toMap(Map.Entry::getKey, e -> e.getValue().length())));

        ProbeQualityProfile probeQualityProfile = new ProbeQualityProfile(
                Map.of(
                        "1", List.of(
                                // Acceptable quality
                                new ProbeQualityWindow(new BaseRegion(1, 10), 1.0f),
                                new ProbeQualityWindow(new BaseRegion(6, 15), 1.0f),
                                new ProbeQualityWindow(new BaseRegion(11, 20), 1.0f),
                                new ProbeQualityWindow(new BaseRegion(16, 25), 1.0f),
                                // Rejected quality
                                new ProbeQualityWindow(new BaseRegion(21, 30), 0.1f),
                                new ProbeQualityWindow(new BaseRegion(26, 35), 0.1f)
                        )
                ),
                10, 5);

        mProbeEvaluator = new ProbeEvaluator(refGenome, probeQualityProfile, null);
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
        Probe probe = new Probe(new ChrBaseRegion("1", 1, 10), METADATA);
        Probe evalProbe = mProbeEvaluator.evaluateProbe(probe, CRITERIA);
        assertNotSame(probe, evalProbe);
        assertEquals(probe.metadata(), evalProbe.metadata());
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
    public void testEvaluateProbesRejectGc()
    {
        Probe probe = new Probe(new ChrBaseRegion("1", 11, 20), METADATA);
        Probe evalProbe = mProbeEvaluator.evaluateProbe(probe, CRITERIA);
        assertNotSame(probe, evalProbe);
        assertEquals(probe.metadata(), evalProbe.metadata());
        assertEquals(CRITERIA, evalProbe.evalCriteria());
        assertEquals("AAAAAAAAAA", evalProbe.sequence());
        assertNotNull(evalProbe.gcContent());
        assertEquals(0, evalProbe.gcContent(), EPSILON);
        assertNotNull(evalProbe.rejectionReason());
        assertTrue(evalProbe.rejected());
        assertFalse(evalProbe.accepted());
        assertTrue(evalProbe.rejectionReason().toLowerCase().contains("gc"));
    }

    @Test
    public void testEvaluateProbesRejectQuality()
    {
        Probe probe = new Probe(new ChrBaseRegion("1", 21, 30), METADATA);
        Probe evalProbe = mProbeEvaluator.evaluateProbe(probe, CRITERIA);
        assertNotSame(probe, evalProbe);
        assertEquals(probe.metadata(), evalProbe.metadata());
        assertEquals(CRITERIA, evalProbe.evalCriteria());
        assertNotNull(evalProbe.qualityScore());
        assertEquals(0.1, evalProbe.qualityScore(), EPSILON);
        assertNotNull(evalProbe.rejectionReason());
        assertTrue(evalProbe.rejected());
        assertFalse(evalProbe.accepted());
        assertTrue(evalProbe.rejectionReason().toLowerCase().contains("qs"));
    }
}
