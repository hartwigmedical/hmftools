package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

public record Probe(
        TargetRegion target,
        ChrBaseRegion region,
        // null if the probe hasn't been evaluated yet.
        @Nullable ProbeEvaluator.Criteria evalCriteria,
        // null if the probe is acceptable.
        @Nullable String rejectionReason,
        // All the following fields may be null if the probe was rejected, and the field didn't need to be checked.
        @Nullable String sequence,
        @Nullable Double qualityScore,
        @Nullable Double gcContent
)
{
    public Probe
    {
        if(!region.chromosome().equals(target.region().chromosome()))
        {
            // Not expecting this to ever occur but should check it to be safe.
            throw new IllegalArgumentException("Probe region and target region should have the same chromosome");
        }
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        if(sequence != null && sequence.length() != region.baseLength())
        {
            throw new IllegalArgumentException("sequence length should match probe region length");
        }
    }

    public Probe(final TargetRegion target, final ChrBaseRegion probeRegion)
    {
        this(target, probeRegion, null, null, null, null, null);
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
        return new Probe(target, region, value, rejectionReason, sequence, qualityScore, gcContent);
    }

    public Probe withRejectionReason(final String value)
    {
        return new Probe(target, region, evalCriteria, value, sequence, qualityScore, gcContent);
    }

    public Probe withSequence(String value)
    {
        return new Probe(target, region, evalCriteria, rejectionReason, value, qualityScore, gcContent);
    }

    public Probe withQualityScore(double value)
    {
        return new Probe(target, region, evalCriteria, rejectionReason, sequence, value, gcContent);
    }

    public Probe withGcContent(double value)
    {
        return new Probe(target, region, evalCriteria, rejectionReason, sequence, qualityScore, value);
    }
}
