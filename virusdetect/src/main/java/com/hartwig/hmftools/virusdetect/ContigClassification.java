package com.hartwig.hmftools.virusdetect;

import java.util.List;

import org.jetbrains.annotations.Nullable;

// One contig's representative-selection verdict.
// Fields fill progressively depending on how far the contig got through the selection process.
public record ContigClassification(
        ViralContig contig,
        ContigFilterStatus filterStatus,
        // Rank by read votes within the oncology group; set for coverage survivors.
        @Nullable Integer votesRank,
        // Read-vote share over all group contigs, and over only the filter-passing candidates. Set for all contigs.
        @Nullable Double preFilterVoteShare,
        @Nullable Double postFilterVoteShare,
        // Set for selection candidates only.
        @Nullable Double voteShareRatio,
        @Nullable Boolean comparable,
        // Ranks of the candidate contigs this one challenges, and is challenged by.
        List<Integer> challengesRanks,
        List<Integer> challengedByRanks,
        @Nullable ContigRole role,
        @Nullable OncologyGroupOutcome oncologyGroupOutcome,
        @Nullable OncologyGroupSubOutcome oncologyGroupSubOutcome)
{
}
