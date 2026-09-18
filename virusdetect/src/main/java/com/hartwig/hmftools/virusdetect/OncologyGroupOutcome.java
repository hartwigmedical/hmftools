package com.hartwig.hmftools.virusdetect;

// The specific outcome of representative contig selection for an oncology group, refining OncologyGroupResolution.
public enum OncologyGroupOutcome
{
    // No candidates: every contig failed the prefilter.
    NO_CANDIDATES(OncologyGroupResolution.NO_CANDIDATES),
    // Resolved: trivially because the oncology group had a single candidate contig.
    ONE_CANDIDATE(OncologyGroupResolution.RESOLVED),
    // Resolved: one clear representative resolved among several candidates.
    RESOLVED_CANDIDATES(OncologyGroupResolution.RESOLVED),
    // Unresolved: a low-abundance contig challenges an abundant one (a possible hidden strain).
    MINOR_RIVAL(OncologyGroupResolution.UNRESOLVED),
    // Unresolved: two abundant contigs challenge each other (a co-infection signature).
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
