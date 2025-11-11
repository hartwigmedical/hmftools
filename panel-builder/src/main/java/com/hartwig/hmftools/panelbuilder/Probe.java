package com.hartwig.hmftools.panelbuilder;

import org.jetbrains.annotations.Nullable;

public record Probe(
        SequenceDefinition definition,
        // null if the probe hasn't been evaluated yet or the probe was rejected by another criteria.
        @Nullable String sequence,
        // Range of bases in the probe sequence which are the intended target of this probe.
        TargetedRange targetedRange,
        TargetMetadata metadata,
        // null if the probe hasn't been evaluated yet.
        @Nullable ProbeEvaluator.Criteria evaluationCriteria,
        // null if the probe hasn't been evaluated yet.
        @Nullable EvaluationResult evaluationResult,
        // null if the probe hasn't been evaluated yet or the probe was rejected by another criteria.
        @Nullable Double qualityScore,
        // null if the probe hasn't been evaluated yet or the probe was rejected by another criteria.
        @Nullable Double gcContent
)
{
    public Probe
    {
        if(sequence != null && definition.baseLength() != sequence.length())
        {
            throw new IllegalArgumentException("sequence length should match definition length");
        }
        if(targetedRange.baseLength() > definition.baseLength())
        {
            throw new IllegalArgumentException("targetedRange must not be larger than the probe");
        }
        if(evaluationResult != null && evaluationCriteria == null)
        {
            throw new IllegalArgumentException("evaluationCriteria must be set if evaluationResult is set");
        }
        if(gcContent != null && !(gcContent >= 0 && gcContent <= 1))
        {
            throw new IllegalArgumentException("gcContent should be between 0 and 1");
        }
        if(qualityScore != null && !(qualityScore >= 0 && qualityScore <= 1))
        {
            throw new IllegalArgumentException("qualityScore should be between 0 and 1");
        }
    }

    public Probe(final SequenceDefinition definition, final TargetedRange targetedRange, final TargetMetadata metadata)
    {
        this(definition, null, targetedRange, metadata, null, null, null, null);
    }

    public boolean evaluated()
    {
        return evaluationResult != null;
    }

    public boolean accepted()
    {
        return evaluated() && evaluationResult.accepted();
    }

    public boolean rejected()
    {
        return evaluated() && evaluationResult.rejected();
    }

    public Probe withSequence(final String value)
    {
        if(sequence != null)
        {
            throw new IllegalArgumentException("sequence already set");
        }
        return new Probe(definition, value, targetedRange, metadata, evaluationCriteria, evaluationResult, qualityScore, gcContent);
    }

    public Probe withEvaluationCriteria(final ProbeEvaluator.Criteria value)
    {
        if(evaluationCriteria != null)
        {
            throw new IllegalArgumentException("evaluationCriteria already set");
        }
        return new Probe(definition, sequence, targetedRange, metadata, value, evaluationResult, qualityScore, gcContent);
    }

    public Probe withEvaluationResult(final EvaluationResult value)
    {
        if(evaluationResult != null)
        {
            throw new IllegalArgumentException("evaluationResult already set");
        }
        return new Probe(definition, sequence, targetedRange, metadata, evaluationCriteria, value, qualityScore, gcContent);
    }

    public Probe withQualityScore(double value)
    {
        if(qualityScore != null)
        {
            throw new IllegalArgumentException("qualityScore already set");
        }
        return new Probe(definition, sequence, targetedRange, metadata, evaluationCriteria, evaluationResult, value, gcContent);
    }

    public Probe withGcContent(double value)
    {
        if(gcContent != null)
        {
            throw new IllegalArgumentException("gcContent already set");
        }
        return new Probe(definition, sequence, targetedRange, metadata, evaluationCriteria, evaluationResult, qualityScore, value);
    }
}
