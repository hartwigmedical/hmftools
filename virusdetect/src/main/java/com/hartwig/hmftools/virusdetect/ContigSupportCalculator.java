package com.hartwig.hmftools.virusdetect;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.virusdetect.VirusConstants.VOTE_CORRECT_BASE_PROBABILITY;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

// Per-contig support statistics over aligned reads: each read counted once per contig, by its single best alignment
// there. Also attributes each read across the contigs it hits as read votes, quantifying strain support within a virus.
public class ContigSupportCalculator
{
    // Held rather than read from the constant directly so tests can inject a controllable value.
    private final double mCorrectBaseProbability;

    public ContigSupportCalculator()
    {
        this(VOTE_CORRECT_BASE_PROBABILITY);
    }

    ContigSupportCalculator(double correctBaseProbability)
    {
        mCorrectBaseProbability = correctBaseProbability;
    }

    public Map<ViralContig, ContigSupport> compute(ViralAlignments viralAlignments)
    {
        Map<ViralContig, ContigAccumulator> accumulators = new HashMap<>();
        viralAlignments.reads().forEach(read -> accumulateRead(read, accumulators));

        // A contig whose every alignment straddled the origin retains no read, but the drop is still worth reporting.
        Map<ViralContig, Integer> originClippedReads = viralAlignments.originClippedReads();
        originClippedReads.keySet().forEach(contig -> accumulators.computeIfAbsent(contig, k -> new ContigAccumulator()));

        // Depth is the costly reduction, so it is taken once here: the group presence decision and each contig's
        // record both need it.
        Map<ViralContig, int[]> depths = accumulators.entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> entry.getValue().calculateDepth(entry.getKey())));

        Set<OncologyGroup> presentGroups = computePresentGroups(depths);
        return accumulators.entrySet().stream().collect(toMap(
                Map.Entry::getKey, entry -> entry.getValue().toContigSupport(
                        entry.getKey(), depths.get(entry.getKey()), originClippedReads.getOrDefault(entry.getKey(), 0),
                        presentGroups.contains(entry.getKey().oncologyGroup()), viralAlignments.meanReadLength())));
    }

    private static Set<OncologyGroup> computePresentGroups(Map<ViralContig, int[]> depths)
    {
        return depths.entrySet().stream()
                .filter(entry -> ContigPrefilter.establishesGroupPresence(coverageFraction(entry.getKey(), entry.getValue())))
                .map(entry -> entry.getKey().oncologyGroup())
                .collect(toSet());
    }

    private static double coverageFraction(ViralContig contig, int[] depth)
    {
        return ContigSupport.coverageFraction(coveredBases(depth), contig);
    }

    private static int coveredBases(int[] depth)
    {
        return (int) Arrays.stream(depth).filter(d -> d > 0).count();
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

    // Accumulates one contig's reads (each read's best alignment there), then reduces them to a ContigSupport.
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

        private ContigSupport toContigSupport(
                ViralContig contig, int[] depth, int originClippedReads, boolean groupPresent, double meanReadLength)
        {
            int coveredBases = coveredBases(depth);
            ContigFilterStatus filterStatus = ContigPrefilter.status(
                    contig, ContigSupport.coverageFraction(coveredBases, contig), mVotes, groupPresent, meanReadLength);
            int multiAlignReads = (int) mAlignmentCounts.stream().filter(count -> count > 1).count();

            // With no retained read there is no distribution to summarise, unlike depth which is zero everywhere.
            boolean hasReads = !mAlignments.isEmpty();
            SummaryStats alignPerRead = hasReads ? SummaryStats.from(mAlignmentCounts) : null;
            SummaryStats alignerScore = hasReads
                    ? SummaryStats.from(mAlignments.stream().mapToInt(ViralAlignment::alignerScore).toArray()) : null;

            return new ContigSupport(
                    contig, filterStatus, mAlignments.size(), multiAlignReads, alignPerRead, originClippedReads,
                    coveredBases, SummaryStats.from(depth), alignerScore, mVotes);
        }

        private int[] calculateDepth(ViralContig contig)
        {
            return ContigSupportCalculator.calculateDepth(contig.length(), mAlignments);
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
