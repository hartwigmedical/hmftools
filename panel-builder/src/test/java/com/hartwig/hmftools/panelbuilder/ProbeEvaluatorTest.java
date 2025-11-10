package com.hartwig.hmftools.panelbuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class ProbeEvaluatorTest
{
    private static final ProbeEvaluator.Criteria CRITERIA = new ProbeEvaluator.Criteria(0.8, 0.5, 0.1);

    // Dummy data, doesn't really matter.
    private static final SequenceDefinition DEFINITION = SequenceDefinition.singleRegion(new ChrBaseRegion("1", 1, 10));
    private static final TargetedRange TARGETED_RANGE = TargetedRange.wholeRegion(DEFINITION.baseLength());
    private static final TargetMetadata METADATA = new TargetMetadata(TargetMetadata.Type.CUSTOM_REGION, "test");

    private static final double EPSILON = 1e-6;

    private static Probe probe(final String sequence, Double qualityScore, Double gcContent)
    {
        return new Probe(DEFINITION, sequence, TARGETED_RANGE, METADATA, null, null, qualityScore, gcContent);
    }

    @Test
    public void testCriteriaGcContentMin()
    {
        assertEquals(0.4, CRITERIA.gcContentMin(), EPSILON);
    }

    @Test
    public void testCriteriaGcContentMax()
    {
        assertEquals(0.6, CRITERIA.gcContentMax(), EPSILON);
    }

    @Test
    public void testEvaluateProbesAcceptable()
    {
        String sequence = "ACGTACGTAC";
        double qualityScore = 1.0;
        double gcContent = 0.5;
        Probe probe = probe(sequence, qualityScore, gcContent);
        Probe evalProbe = ProbeEvaluator.evaluateProbe(probe, CRITERIA);
        assertNotSame(probe, evalProbe);
        assertEquals(probe.metadata(), evalProbe.metadata());
        assertEquals(CRITERIA, evalProbe.evalCriteria());
        assertEquals(qualityScore, evalProbe.qualityScore(), EPSILON);
        assertEquals(sequence, evalProbe.sequence());
        assertEquals(gcContent, evalProbe.gcContent(), EPSILON);
        assertNull(evalProbe.rejectionReason());
        assertTrue(evalProbe.accepted());
        assertFalse(evalProbe.rejected());
    }

    @Test
    public void testEvaluateProbesRejectSequence()
    {
        String sequence = "AAANNNAAAA";
        Probe probe = probe(sequence, 1.0, 0.5);
        Probe evalProbe = ProbeEvaluator.evaluateProbe(probe, CRITERIA);
        assertNotSame(probe, evalProbe);
        assertEquals(probe.metadata(), evalProbe.metadata());
        assertEquals(CRITERIA, evalProbe.evalCriteria());
        assertEquals(sequence, evalProbe.sequence());
        assertNotNull(evalProbe.rejectionReason());
        assertTrue(evalProbe.rejected());
        assertFalse(evalProbe.accepted());
        assertTrue(evalProbe.rejectionReason().toLowerCase().contains("seq"));
    }

    @Test
    public void testEvaluateProbesRejectGc()
    {
        Probe probe = probe("AAAAAAAAAA", 1.0, 0.0);
        Probe evalProbe = ProbeEvaluator.evaluateProbe(probe, CRITERIA);
        assertNotSame(probe, evalProbe);
        assertEquals(probe.metadata(), evalProbe.metadata());
        assertEquals(CRITERIA, evalProbe.evalCriteria());
        assertEquals(0, evalProbe.gcContent(), EPSILON);
        assertNotNull(evalProbe.rejectionReason());
        assertTrue(evalProbe.rejected());
        assertFalse(evalProbe.accepted());
        assertTrue(evalProbe.rejectionReason().toLowerCase().contains("gc"));
    }

    @Test
    public void testEvaluateProbesRejectQuality()
    {
        Probe probe = probe("ACGTACGTAC", 0.1, 0.5);
        Probe evalProbe = ProbeEvaluator.evaluateProbe(probe, CRITERIA);
        assertNotSame(probe, evalProbe);
        assertEquals(probe.metadata(), evalProbe.metadata());
        assertEquals(CRITERIA, evalProbe.evalCriteria());
        assertEquals(0.1, evalProbe.qualityScore(), EPSILON);
        assertNotNull(evalProbe.rejectionReason());
        assertTrue(evalProbe.rejected());
        assertFalse(evalProbe.accepted());
        assertTrue(evalProbe.rejectionReason().toLowerCase().contains("qs"));
    }

    // TODO: test annotating attributes
}
