package com.hartwig.hmftools.panelbuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.function.Function;
import java.util.stream.Stream;

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

    // Attributes are pre-set in test cases, so no need for these functions to modify the probe. Simplifies the test setup.
    private final ProbeEvaluator mEvaluator = new ProbeEvaluator(Function.identity(), Function.identity(), Function.identity());

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
        Probe probe = probe(sequence, qualityScore, gcContent).withEvaluationCriteria(CRITERIA);
        Probe evalProbe = mEvaluator.evaluateProbes(Stream.of(probe)).findFirst().orElseThrow();
        assertEquals(probe.metadata(), evalProbe.metadata());
        assertEquals(CRITERIA, evalProbe.evaluationCriteria());
        assertEquals(qualityScore, evalProbe.qualityScore(), EPSILON);
        assertEquals(sequence, evalProbe.sequence());
        assertEquals(gcContent, evalProbe.gcContent(), EPSILON);
        assertEquals(EvaluationResult.accept(), evalProbe.evaluationResult());
        assertTrue(evalProbe.accepted());
        assertFalse(evalProbe.rejected());
    }

    @Test
    public void testEvaluateProbesRejectSequence()
    {
        String sequence = "AAANNNAAAA";
        Probe probe = probe(sequence, 1.0, 0.5).withEvaluationCriteria(CRITERIA);
        Probe evalProbe = mEvaluator.evaluateProbes(Stream.of(probe)).findFirst().orElseThrow();
        assertEquals(probe.metadata(), evalProbe.metadata());
        assertEquals(CRITERIA, evalProbe.evaluationCriteria());
        assertEquals(sequence, evalProbe.sequence());
        assertTrue(evalProbe.rejected());
        assertFalse(evalProbe.accepted());
        assertEquals(EvaluationResult.reject("sequence"), evalProbe.evaluationResult());
    }

    @Test
    public void testEvaluateProbesRejectGc()
    {
        double gcContent = 0.0;
        Probe probe = probe("AAAAAAAAAA", 1.0, gcContent).withEvaluationCriteria(CRITERIA);
        Probe evalProbe = mEvaluator.evaluateProbes(Stream.of(probe)).findFirst().orElseThrow();
        assertEquals(probe.metadata(), evalProbe.metadata());
        assertEquals(CRITERIA, evalProbe.evaluationCriteria());
        assertEquals(gcContent, evalProbe.gcContent(), EPSILON);
        assertTrue(evalProbe.rejected());
        assertFalse(evalProbe.accepted());
        assertEquals(EvaluationResult.reject("GC"), evalProbe.evaluationResult());
    }

    @Test
    public void testEvaluateProbesRejectQuality()
    {
        double qualityScore = 0.1;
        Probe probe = probe("ACGTACGTAC", qualityScore, 0.5).withEvaluationCriteria(CRITERIA);
        Probe evalProbe = mEvaluator.evaluateProbes(Stream.of(probe)).findFirst().orElseThrow();
        assertEquals(probe.metadata(), evalProbe.metadata());
        assertEquals(CRITERIA, evalProbe.evaluationCriteria());
        assertEquals(qualityScore, evalProbe.qualityScore(), EPSILON);
        assertTrue(evalProbe.rejected());
        assertFalse(evalProbe.accepted());
        assertEquals(EvaluationResult.reject("QS"), evalProbe.evaluationResult());
    }

    // TODO: test annotating attributes
}
