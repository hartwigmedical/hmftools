package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.Utils.isDnaSequenceNormal;

import org.jetbrains.annotations.Nullable;

public record Probe(
        SequenceDefinition definition,
        String sequence,
        TargetMetadata metadata,
        // null if the probe hasn't been evaluated yet.
        @Nullable ProbeEvaluator.Criteria evalCriteria,
        // null if the probe is acceptable.
        @Nullable String rejectionReason,
        double qualityScore,
        double gcContent
)
{
    public Probe
    {
        if(definition.baseLength() != sequence.length())
        {
            throw new IllegalArgumentException("sequence length should match definition length");
        }
        if(!isDnaSequenceNormal(sequence))
        {
            throw new IllegalArgumentException("sequence must only contain normal DNA bases");
        }
        if(rejectionReason != null && (rejectionReason.isBlank()))
        {
            throw new IllegalArgumentException("rejectionReason should not be blank");
        }
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

    public Probe withEvalCriteria(final ProbeEvaluator.Criteria value)
    {
        return new Probe(definition, sequence, metadata, value, rejectionReason, qualityScore, gcContent);
    }

    public Probe withRejectionReason(final String value)
    {
        return new Probe(definition, sequence, metadata, evalCriteria, value, qualityScore, gcContent);
    }
}
