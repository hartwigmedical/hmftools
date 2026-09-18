package com.hartwig.hmftools.virusdetect;

import java.util.List;

// Result of a candidate contig in representative contig selection.
public record CandidateSelectionResult(
        // Rank by read votes among the candidates. 1 is the highest.
        int votesRank,
        // Indicates if read votes are near the group's highest - abundant enough to be a rival.
        boolean comparable,
        // Votes-ranks of the candidates this one challenges, and of those which challenge it.
        List<Integer> challengesRanks,
        List<Integer> challengedByRanks,
        ContigRole role
)
{
}
