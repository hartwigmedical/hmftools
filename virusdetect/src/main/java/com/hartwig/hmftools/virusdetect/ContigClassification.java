package com.hartwig.hmftools.virusdetect;

import java.util.List;

import org.jetbrains.annotations.Nullable;

// One contig's representative-selection verdict.
// Fields fill progressively depending on how far the contig got through the selection process.
public record ContigClassification(
        String contig,
        String oncologyGroup,
        ContigFilterStatus filterStatus,
        // Set for coverage survivors: rank by read votes within the oncology group, and share of its total votes.
        @Nullable Integer votesRank,
        @Nullable Double voteShare,
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
