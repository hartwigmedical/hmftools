package com.hartwig.hmftools.virusdetect;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.virusdetect.VirusConstants.VOTE_CORRECT_BASE_PROBABILITY;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// Per-contig support statistics over aligned reads: each read counted once per contig, by its single best alignment
// there. Also attributes each read across the contigs it hits as read votes, quantifying strain support within a virus.
public class ContigStatsCalculator
{
    // Held rather than read from the constant directly so tests can inject a controllable value.
    private final double mCorrectBaseProbability;

    public ContigStatsCalculator()
    {
        this(VOTE_CORRECT_BASE_PROBABILITY);
    }

    ContigStatsCalculator(double correctBaseProbability)
    {
        mCorrectBaseProbability = correctBaseProbability;
    }

    public Map<ViralContig, ContigStats> compute(ViralAlignments viralAlignments)
    {
        Map<ViralContig, ContigAccumulator> accumulators = new HashMap<>();
        viralAlignments.reads().forEach(read -> accumulateRead(read, accumulators));

        // A contig whose every alignment straddled the origin retains no read, but the drop is still worth reporting.
        Map<ViralContig, Integer> originClippedReads = viralAlignments.originClippedReads();
        originClippedReads.keySet().forEach(contig -> accumulators.computeIfAbsent(contig, k -> new ContigAccumulator()));

        return accumulators.entrySet().stream().collect(toMap(
                Map.Entry::getKey, entry -> entry.getValue().toContigStats(
                        entry.getKey(), originClippedReads.getOrDefault(entry.getKey(), 0))));
    }

    // Folds one read into the per-contig accumulators: its best alignment (and alignment count) on each contig it hits,
    // plus its cross-contig vote split.
    private void accumulateRead(ReadAlignments read, Map<ViralContig, ContigAccumulator> accumulators)
    {
        read.hits().forEach((contig, hit) ->
                accumulators.computeIfAbsent(contig, k -> new ContigAccumulator()).addRead(hit.best(), hit.alignmentCount()));

        addVotes(read, accumulators);
    }

    // A read's vote splits across the contigs it hits by how well each explains it: every extra base a contig fails to
    // explain multiplies its share by the (pessimistic) chance a base is right, so the closest contig wins most.
    private void addVotes(ReadAlignments read, Map<ViralContig, ContigAccumulator> accumulators)
    {
        int minDivergence = read.minDivergence();
        Map<ViralContig, Double> contigWeights = read.hits().entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> voteWeight(entry.getValue().divergence(), minDivergence)));

        double totalWeight = contigWeights.values().stream().mapToDouble(Double::doubleValue).sum();
        contigWeights.forEach((contig, weight) -> accumulators.get(contig).addVote(weight / totalWeight));
    }

    private double voteWeight(int divergence, int minDivergence)
    {
        return Math.pow(mCorrectBaseProbability, divergence - minDivergence);
    }

    // Accumulates one contig's reads (each read's best alignment there), then reduces them to a ContigStats.
    private static class ContigAccumulator
    {
        private final List<ViralAlignment> mAlignments = new ArrayList<>();
        private final List<Integer> mAlignmentCounts = new ArrayList<>();
        private double mVotes;

        private void addRead(ViralAlignment best, int alignmentCount)
        {
            mAlignments.add(best);
            mAlignmentCounts.add(alignmentCount);
        }

        private void addVote(double vote)
        {
            mVotes += vote;
        }

        private ContigStats toContigStats(ViralContig contig, int originClippedReads)
        {
            int[] depth = calculateDepth(contig.length(), mAlignments);
            int coveredBases = (int) Arrays.stream(depth).filter(d -> d > 0).count();
            int multiAlignReads = (int) mAlignmentCounts.stream().filter(count -> count > 1).count();

            // With no retained read there is no distribution to summarise, unlike depth which is zero everywhere.
            boolean hasReads = !mAlignments.isEmpty();
            SummaryStats alignPerRead = hasReads ? SummaryStats.from(mAlignmentCounts) : null;
            SummaryStats alignerScore = hasReads
                    ? SummaryStats.from(mAlignments.stream().mapToInt(ViralAlignment::alignerScore).toArray()) : null;

            return new ContigStats(
                    contig, mAlignments.size(), multiAlignReads, alignPerRead, originClippedReads,
                    coveredBases, SummaryStats.from(depth), alignerScore, mVotes);
        }
    }

    private static int[] calculateDepth(int contigLength, Collection<ViralAlignment> alignments)
    {
        int[] depth = new int[contigLength];
        for(ViralAlignment alignment : alignments)
        {
            for(AlignedInterval block : alignment.alignedIntervals())
            {
                int start = max(0, block.referenceStart() - 1);   // intervals are 1-based
                int end = min(contigLength, block.referenceStart() - 1 + block.length());
                for(int position = start; position < end; ++position)
                {
                    ++depth[position];
                }
            }
        }
        return depth;
    }
}
