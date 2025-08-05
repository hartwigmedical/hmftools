package com.hartwig.hmftools.panelbuilder;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

public record Probe(
        // null if the probe sequence doesn't correspond to a single location in the reference genome.
        // This is the case when designing probes to capture variants.
        @Nullable ChrBaseRegion region,
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
        if(region != null && !region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        if(region != null && sequence.length() != region.baseLength())
        {
            throw new IllegalArgumentException("sequence length should match probe region length");
        }
        if (!sequence.matches("[ACGT]*")) {
            throw new IllegalArgumentException("sequence must only contain valid DNA bases");
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
        return new Probe(region, sequence, metadata, value, rejectionReason, qualityScore, gcContent);
    }

    public Probe withRejectionReason(final String value)
    {
        return new Probe(region, sequence, metadata, evalCriteria, value, qualityScore, gcContent);
    }
}
