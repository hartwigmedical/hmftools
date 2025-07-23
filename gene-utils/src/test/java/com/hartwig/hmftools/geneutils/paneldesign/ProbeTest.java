package com.hartwig.hmftools.geneutils.paneldesign;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class ProbeTest
{
    private static final TargetRegion TARGET_REGION = new TargetRegion(
            new ChrBaseRegion("1", 100, 200),
            new TargetMetadata(TargetMetadata.Type.CUSTOM, "test"));
    private static final ChrBaseRegion PROBE_REGION = new ChrBaseRegion("1", 100, 109);
    private static final String PROBE_SEQUENCE = "ACGTACGTAC";
    private static final ProbeEvaluator.Criteria EVAL_CRITERIA = new ProbeEvaluator.Criteria(0.8, 0.45, 0.1);

    @Test
    public void testConstructor()
    {
        Probe probe = new Probe(TARGET_REGION, PROBE_REGION);
        assertEquals(
                new Probe(TARGET_REGION, PROBE_REGION, null, null, null, null, null),
                probe);
        assertEquals(TARGET_REGION, probe.target());
        assertEquals(PROBE_REGION, probe.region());
        assertNull(probe.evalCriteria());
        assertNull(probe.sequence());
        assertNull(probe.qualityScore());
        assertNull(probe.gcContent());
    }

    @Test
    public void testWith()
    {
        Probe probe = new Probe(TARGET_REGION, PROBE_REGION);

        probe = probe.withEvalCriteria(EVAL_CRITERIA);
        assertEquals(EVAL_CRITERIA, probe.evalCriteria());

        String rejectionReason = "rejected";
        probe = probe.withRejectionReason(rejectionReason);
        assertEquals(rejectionReason, probe.rejectionReason());

        probe = probe.withSequence(PROBE_SEQUENCE);
        assertEquals(PROBE_SEQUENCE, probe.sequence());

        double qualityScore = 0.1;
        probe = probe.withQualityScore(qualityScore);
        assertEquals(qualityScore, probe.qualityScore().doubleValue(), 0);

        double gcContent = 0.2;
        probe = probe.withGcContent(gcContent);
        assertEquals(gcContent, probe.gcContent().doubleValue(), 0);

        assertEquals(new Probe(TARGET_REGION, PROBE_REGION, EVAL_CRITERIA, rejectionReason, PROBE_SEQUENCE, qualityScore, gcContent), probe);
    }

    @Test
    public void testEvalAccepted()
    {
        Probe probe = new Probe(TARGET_REGION, PROBE_REGION);
        assertFalse(probe.evaluated());
        assertFalse(probe.accepted());
        assertFalse(probe.rejected());

        probe = probe.withEvalCriteria(EVAL_CRITERIA);
        probe = probe.withQualityScore(1.0);
        probe = probe.withGcContent(0.45);
        probe = probe.withSequence(PROBE_SEQUENCE);
        assertTrue(probe.evaluated());
        assertTrue(probe.accepted());
        assertFalse(probe.rejected());
    }

    @Test
    public void testEvalRejected()
    {
        Probe probe = new Probe(TARGET_REGION, PROBE_REGION);
        assertFalse(probe.evaluated());
        assertFalse(probe.accepted());
        assertFalse(probe.rejected());

        probe = probe.withEvalCriteria(EVAL_CRITERIA);
        probe = probe.withQualityScore(0.1);
        probe = probe.withRejectionReason("quality");
        assertTrue(probe.evaluated());
        assertFalse(probe.accepted());
        assertTrue(probe.rejected());
    }
}
