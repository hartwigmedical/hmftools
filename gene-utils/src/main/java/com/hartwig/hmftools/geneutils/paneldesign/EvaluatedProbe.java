package com.hartwig.hmftools.geneutils.paneldesign;

public record EvaluatedProbe(
        CandidateProbe candidate,
        // null if the probe is acceptable.
        String rejectionReason,
        // null if not checked.
        Double qualityScore,
        // null if not checked.
        Double gcContent
)
{
    public EvaluatedProbe(final CandidateProbe probe)
    {
        this(probe, null, null, null);
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
        return new EvaluatedProbe(candidate, value, qualityScore, gcContent);
    }

    public EvaluatedProbe withQualityScore(double value)
    {
        return new EvaluatedProbe(candidate, rejectionReason, value, gcContent);
    }

    public EvaluatedProbe withGcContent(double value)
    {
        return new EvaluatedProbe(candidate, rejectionReason, qualityScore, value);
    }
}
