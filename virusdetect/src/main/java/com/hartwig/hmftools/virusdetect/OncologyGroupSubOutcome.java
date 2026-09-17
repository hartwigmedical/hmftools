package com.hartwig.hmftools.virusdetect;

// The specific reason an oncology group resolved or not, refining OncologyGroupOutcome.
public enum OncologyGroupSubOutcome
{
    // Resolved: trivially because the oncology group had a single candidate contig.
    ONE_CANDIDATE(OncologyGroupOutcome.RESOLVED),
    // Resolved: one clear representative resolved among several candidates.
    RESOLVED_CANDIDATES(OncologyGroupOutcome.RESOLVED),
    // Unresolved: a low-abundance contig challenges an abundant one (a possible hidden strain).
    MINOR_RIVAL(OncologyGroupOutcome.UNRESOLVED),
    // Unresolved: two abundant contigs challenge each other (a co-infection signature).
    MUTUAL(OncologyGroupOutcome.UNRESOLVED),
    // Unresolved: the challenge relation loops, with no clear top.
    CYCLE(OncologyGroupOutcome.UNRESOLVED);

    private final OncologyGroupOutcome mOutcome;

    private OncologyGroupSubOutcome(OncologyGroupOutcome outcome)
    {
        mOutcome = outcome;
    }

    public OncologyGroupOutcome outcome()
    {
        return mOutcome;
    }
}
