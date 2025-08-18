package com.hartwig.hmftools.panelbuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class ProbeTest
{
    private static final TargetMetadata METADATA = new TargetMetadata(TargetMetadata.Type.CUSTOM, "test");
    private static final ProbeTarget TARGET = ProbeTarget.exactRegion(new ChrBaseRegion("1", 100, 129));
    private static final String SEQUENCE = "ACGTACGTACACGTACGTACACGTACGTAC";
    private static final ProbeEvaluator.Criteria EVAL_CRITERIA = new ProbeEvaluator.Criteria(0.8, 0.45, 0.1);

    @Test
    public void testWith()
    {
        double qualityScore = 0.1;
        double gcContent = 0.2;
        Probe probe = new Probe(TARGET, SEQUENCE, METADATA, null, null, qualityScore, gcContent);

        probe = probe.withEvalCriteria(EVAL_CRITERIA);
        assertEquals(EVAL_CRITERIA, probe.evalCriteria());

        String rejectionReason = "rejected";
        probe = probe.withRejectionReason(rejectionReason);
        assertEquals(rejectionReason, probe.rejectionReason());

        assertEquals(new Probe(TARGET, SEQUENCE, METADATA, EVAL_CRITERIA, rejectionReason, qualityScore, gcContent), probe);
    }

    @Test
    public void testEvalAccepted()
    {
        Probe probe = new Probe(TARGET, SEQUENCE, METADATA, null, null, 1.0, 0.45);

        assertFalse(probe.evaluated());
        assertFalse(probe.accepted());
        assertFalse(probe.rejected());

        probe = probe.withEvalCriteria(EVAL_CRITERIA);
        assertTrue(probe.evaluated());
        assertTrue(probe.accepted());
        assertFalse(probe.rejected());
    }

    @Test
    public void testEvalRejected()
    {
        Probe probe = new Probe(TARGET, SEQUENCE, METADATA, null, null, 0.1, 0.45);

        assertFalse(probe.evaluated());
        assertFalse(probe.accepted());
        assertFalse(probe.rejected());

        probe = probe.withEvalCriteria(EVAL_CRITERIA);
        probe = probe.withRejectionReason("QS");
        assertTrue(probe.evaluated());
        assertFalse(probe.accepted());
        assertTrue(probe.rejected());
    }
}
