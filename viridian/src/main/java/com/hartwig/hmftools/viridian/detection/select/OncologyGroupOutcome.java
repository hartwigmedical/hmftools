package com.hartwig.hmftools.viridian.detection.select;

// The specific outcome of representative genome selection for an oncology group, refining OncologyGroupResolution.
public enum OncologyGroupOutcome
{
    // No candidates: every genome failed the support filters.
    NO_CANDIDATES(OncologyGroupResolution.NO_CANDIDATES),
    // Resolved: trivially because the oncology group had a single candidate genome.
    ONE_CANDIDATE(OncologyGroupResolution.RESOLVED),
    // Resolved: one clear representative resolved among several candidates.
    RESOLVED_CANDIDATES(OncologyGroupResolution.RESOLVED),
    // Unresolved: a low-abundance genome challenges an abundant one (a possible hidden strain).
    MINOR_RIVAL(OncologyGroupResolution.UNRESOLVED),
    // Unresolved: two abundant genomes challenge each other (a co-infection signature).
    MUTUAL(OncologyGroupResolution.UNRESOLVED),
    // Unresolved: the challenge relation loops, with no clear top.
    CYCLE(OncologyGroupResolution.UNRESOLVED);

    private final OncologyGroupResolution mResolution;

    private OncologyGroupOutcome(OncologyGroupResolution resolution)
    {
        mResolution = resolution;
    }

    public OncologyGroupResolution resolution()
    {
        return mResolution;
    }
}
