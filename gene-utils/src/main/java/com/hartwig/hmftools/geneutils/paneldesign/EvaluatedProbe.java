package com.hartwig.hmftools.geneutils.paneldesign;

import org.jetbrains.annotations.Nullable;

public record EvaluatedProbe(
        CandidateProbe candidate,
        ProbeEvalCriteria criteria,
        // null if the probe is acceptable.
        @Nullable String rejectionReason,
        // All the following fields may be null if the probe was rejected, and the field didn't need to be checked.
        @Nullable String sequence,
        @Nullable Double qualityScore,
        @Nullable Double gcContent
)
{
    public EvaluatedProbe(final CandidateProbe probe, final ProbeEvalCriteria criteria)
    {
        this(probe, criteria, null, null, null, null);
    }

    public boolean accepted()
    {
        return rejectionReason == null;
    }

    public boolean rejected()
    {
        return !accepted();
    }

    public EvaluatedProbe withRejectionReason(final String value)
    {
        return new EvaluatedProbe(candidate, criteria, value, sequence, qualityScore, gcContent);
    }

    public EvaluatedProbe withSequence(String value)
    {
        return new EvaluatedProbe(candidate, criteria, rejectionReason, value, qualityScore, gcContent);
    }

    public EvaluatedProbe withQualityScore(double value)
    {
        return new EvaluatedProbe(candidate, criteria, rejectionReason, sequence, value, gcContent);
    }

    public EvaluatedProbe withGcContent(double value)
    {
        return new EvaluatedProbe(candidate, criteria, rejectionReason, sequence, qualityScore, value);
    }
}
