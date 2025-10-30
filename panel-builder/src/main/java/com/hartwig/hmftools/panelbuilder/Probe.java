package com.hartwig.hmftools.panelbuilder;

import org.jetbrains.annotations.Nullable;

public record Probe(
        SequenceDefinition definition,
        // null if the probe hasn't been evaluated yet.
        @Nullable String sequence,
        TargetMetadata metadata,
        // null if the probe hasn't been evaluated yet.
        @Nullable ProbeEvaluator.Criteria evalCriteria,
        // null if the probe is acceptable.
        @Nullable String rejectionReason,
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
        if(rejectionReason != null && rejectionReason.isBlank())
        {
            throw new IllegalArgumentException("rejectionReason should not be blank");
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

    public Probe(final SequenceDefinition definition, final TargetMetadata metadata)
    {
        this(definition, null, metadata, null, null, null, null);
    }

    public boolean evaluated()
    {
        return evalCriteria != null;
    }

    public boolean accepted()
    {
        return evaluated() && rejectionReason == null;
    }

    public boolean rejected()
    {
        return evaluated() && rejectionReason != null;
    }

    public Probe withSequence(final String value)
    {
        if(sequence != null)
        {
            throw new IllegalArgumentException("sequence already set");
        }
        return new Probe(definition, value, metadata, evalCriteria, rejectionReason, qualityScore, gcContent);
    }

    public Probe withEvalCriteria(final ProbeEvaluator.Criteria value)
    {
        if(evalCriteria != null)
        {
            throw new IllegalArgumentException("evalCriteria already set");
        }
        return new Probe(definition, sequence, metadata, value, rejectionReason, qualityScore, gcContent);
    }

    public Probe withRejectionReason(final String value)
    {
        if(rejectionReason != null)
        {
            throw new IllegalArgumentException("rejectionReason already set");
        }
        return new Probe(definition, sequence, metadata, evalCriteria, value, qualityScore, gcContent);
    }

    public Probe withQualityScore(double value)
    {
        if(qualityScore != null)
        {
            throw new IllegalArgumentException("qualityScore already set");
        }
        return new Probe(definition, sequence, metadata, evalCriteria, rejectionReason, value, gcContent);
    }

    public Probe withGcContent(double value)
    {
        if(gcContent != null)
        {
            throw new IllegalArgumentException("gcContent already set");
        }
        return new Probe(definition, sequence, metadata, evalCriteria, rejectionReason, qualityScore, value);
    }
}
