package com.hartwig.hmftools.geneutils.paneldesign;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
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

    private final Probe mProbe = new Probe(TARGET_REGION, PROBE_REGION);

    @Test
    public void testConstructor()
    {
        assertEquals(
                new Probe(TARGET_REGION, PROBE_REGION, null, null, null, null, null),
                mProbe);
        assertEquals(TARGET_REGION, mProbe.target());
        assertEquals(PROBE_REGION, mProbe.region());
        assertNull(mProbe.evalCriteria());
        assertNull(mProbe.sequence());
        assertNull(mProbe.qualityScore());
        assertNull(mProbe.gcContent());
    }

    @Test
    public void testWith()
    {
        Probe probe = mProbe;

        probe = probe.withEvalCriteria(EVAL_CRITERIA);
        assertEquals(EVAL_CRITERIA, probe.evalCriteria());

        String rejectionReason = "rejected";
        probe = probe.withRejectionReason(rejectionReason);
        assertEquals(rejectionReason, probe.rejectionReason());

        probe = probe.withSequence(PROBE_SEQUENCE);
        assertEquals(PROBE_SEQUENCE, probe.sequence());

        double qualityScore = 0.1;
        probe = probe.withQualityScore(qualityScore);
        assertNotNull(probe.qualityScore());
        assertEquals(qualityScore, probe.qualityScore(), 0);

        double gcContent = 0.2;
        probe = probe.withGcContent(gcContent);
        assertNotNull(probe.gcContent());
        assertEquals(gcContent, probe.gcContent(), 0);

        assertEquals(new Probe(TARGET_REGION, PROBE_REGION, EVAL_CRITERIA, rejectionReason, PROBE_SEQUENCE, qualityScore, gcContent), probe);
    }

    @Test
    public void testEvalAccepted()
    {
        Probe probe = mProbe;

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
        Probe probe = mProbe;

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
