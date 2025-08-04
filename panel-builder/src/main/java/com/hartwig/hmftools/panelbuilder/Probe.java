package com.hartwig.hmftools.panelbuilder;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

public record Probe(
        // null if the probe sequence doesn't correspond to a single location in the reference genome.
        // This is the case when designing probes to capture variants.
        @Nullable ChrBaseRegion region,
        // null if region is not null and the probe got rejected before the sequence needed to be checked.
        @Nullable String sequence,
        TargetMetadata metadata,
        // null if the probe hasn't been evaluated yet.
        @Nullable ProbeEvaluator.Criteria evalCriteria,
        // null if the probe is acceptable.
        @Nullable String rejectionReason,
        // All the following fields may be null if the probe was rejected, and the field didn't need to be checked.
        @Nullable Double qualityScore,
        @Nullable Double gcContent
)
{
    public Probe
    {
        if(region == null && sequence == null)
        {
            throw new IllegalArgumentException("Either region or sequence must be specified");
        }
        if(region != null && !region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        if(sequence != null && region != null && sequence.length() != region.baseLength())
        {
            throw new IllegalArgumentException("sequence length should match probe region length");
        }
        if(rejectionReason != null && (rejectionReason.isBlank()))
        {
            throw new IllegalArgumentException("rejectionReason should not be blank");
        }
    }

    public Probe(final ChrBaseRegion region, final TargetMetadata metadata)
    {
        this(region, null, metadata, null, null, null, null);
    }

    public Probe(final String sequence, final TargetMetadata metadata)
    {
        this(null, sequence, metadata, null, null, null, null);
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

    public Probe withSequence(String value)
    {
        return new Probe(region, value, metadata, evalCriteria, rejectionReason, qualityScore, gcContent);
    }

    public Probe withQualityScore(double value)
    {
        return new Probe(region, sequence, metadata, evalCriteria, rejectionReason, value, gcContent);
    }

    public Probe withGcContent(double value)
    {
        return new Probe(region, sequence, metadata, evalCriteria, rejectionReason, qualityScore, value);
    }
}
