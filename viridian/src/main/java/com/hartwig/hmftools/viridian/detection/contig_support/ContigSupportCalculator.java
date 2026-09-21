package com.hartwig.hmftools.viridian.detection.contig_support;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.READ_VOTE_CORRECT_BASE_PROBABILITY;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.common.SummaryStats;
import com.hartwig.hmftools.viridian.detection.read_align.AlignedInterval;
import com.hartwig.hmftools.viridian.detection.read_align.AlignedRead;
import com.hartwig.hmftools.viridian.detection.read_align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.detection.read_align.ViralReadAlignments;
import com.hartwig.hmftools.viridian.reference.ViralContig;

// Per-contig support over aligned reads: each read counted once per contig, by its single best alignment there.
// Also attributes each read across the contigs it hits as read votes, quantifying strain support within a virus,
// and judges from coverage and votes whether the contig is present in the sample at all.
public class ContigSupportCalculator
{
    // Held rather than read from the constant directly so tests can inject a controllable value.
    private final double mCorrectBaseProbability;

    public ContigSupportCalculator(double correctBaseProbability)
    {
        mCorrectBaseProbability = correctBaseProbability;
    }

    public ContigSupportCalculator()
    {
        this(READ_VOTE_CORRECT_BASE_PROBABILITY);
    }

    public List<ContigSupport> compute(ViralReadAlignments viralAlignments)
    {
        Map<ViralContig, ContigAccumulator> accumulators = new HashMap<>();
        viralAlignments.reads().forEach(read -> processRead(read, accumulators));

        // A contig whose every alignment straddled the origin retains no read, but the drop is still worth reporting.
        Map<ViralContig, Integer> originClippedReads = viralAlignments.originClippedReads();
        originClippedReads.keySet().forEach(contig -> accumulators.computeIfAbsent(contig, ContigAccumulator::new));

        return createContigSupports(accumulators, originClippedReads, viralAlignments.meanReadLength());
    }

    private static class ContigAccumulator
    {
        public final ViralContig Contig;
        public final List<ViralReadAlignment> Alignments = new ArrayList<>();
        public final List<Integer> AlignmentCounts = new ArrayList<>();
        public double Votes;

        private ContigAccumulator(ViralContig contig)
        {
            Contig = contig;
        }
    }

    private void processRead(AlignedRead read, Map<ViralContig, ContigAccumulator> accumulators)
    {
        read.hits().forEach((contig, hit) ->
        {
            ContigAccumulator accumulator = accumulators.computeIfAbsent(contig, ContigAccumulator::new);
            accumulator.Alignments.add(hit.best());
            accumulator.AlignmentCounts.add(hit.alignmentCount());
        });

        addVotes(read, accumulators);
    }

    // A read's vote splits across the contigs it hits by how well each explains it: every extra base a contig fails to
    // explain multiplies its share by the (pessimistic) chance a base is right, so the closest contig wins most.
    private void addVotes(AlignedRead read, Map<ViralContig, ContigAccumulator> accumulators)
    {
        int minDivergence = read.minDivergence();
        Map<ViralContig, Double> contigWeights = read.hits().entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> voteWeight(entry.getValue().divergence(), minDivergence)));

        double totalWeight = contigWeights.values().stream().mapToDouble(Double::doubleValue).sum();
        contigWeights.forEach((contig, weight) -> accumulators.get(contig).Votes += weight / totalWeight);
    }

    private double voteWeight(int divergence, int minDivergence)
    {
        return Math.pow(mCorrectBaseProbability, divergence - minDivergence);
    }

    private static List<ContigSupport> createContigSupports(
            Map<ViralContig, ContigAccumulator> accumulators, Map<ViralContig, Integer> originClippedReads,
            double meanReadLength)
    {
        return accumulators.values().stream()
                .collect(groupingBy(accumulator -> accumulator.Contig.oncologyGroup()))
                .values().stream()
                .flatMap(group -> createGroupContigSupports(group, originClippedReads, meanReadLength).stream())
                .toList();
    }

    private static List<ContigSupport> createGroupContigSupports(
            List<ContigAccumulator> group, Map<ViralContig, Integer> originClippedReads, double meanReadLength)
    {
        Map<ViralContig, int[]> depths = group.stream().collect(toMap(
                accumulator -> accumulator.Contig,
                accumulator -> calculateDepth(accumulator.Contig.length(), accumulator.Alignments)));

        Map<ViralContig, Integer> coveredBases = depths.entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> countCoveredBases(entry.getValue())));
        Map<ViralContig, Double> readVotes = group.stream().collect(toMap(
                accumulator -> accumulator.Contig, accumulator -> accumulator.Votes));

        Map<ViralContig, ContigFilterStatus> filterStatuses = ContigSupportFilter.statuses(coveredBases, readVotes, meanReadLength);

        return group.stream().map(accumulator ->
        {
            ViralContig contig = accumulator.Contig;
            return createContigSupport(
                    accumulator, filterStatuses.get(contig), coveredBases.get(contig), SummaryStats.from(depths.get(contig)),
                    originClippedReads.getOrDefault(contig, 0));
        }).toList();
    }

    private static ContigSupport createContigSupport(
            ContigAccumulator accumulator, ContigFilterStatus filterStatus, int coveredBases, SummaryStats depth,
            int originClippedReads)
    {
        int multiAlignReads = (int) accumulator.AlignmentCounts.stream().filter(count -> count > 1).count();

        // Without retained reads there is no distribution to summarise, unlike depth which is zero everywhere.
        boolean hasReads = !accumulator.Alignments.isEmpty();
        SummaryStats alignmentsPerRead = hasReads ? SummaryStats.from(accumulator.AlignmentCounts) : null;
        SummaryStats alignerScore = hasReads
                ? SummaryStats.from(accumulator.Alignments.stream().mapToInt(ViralReadAlignment::alignerScore).toArray()) : null;

        return new ContigSupport(
                accumulator.Contig, filterStatus, accumulator.Alignments.size(), multiAlignReads, alignmentsPerRead,
                originClippedReads, coveredBases, depth, alignerScore, accumulator.Votes);
    }

    private static int countCoveredBases(int[] depth)
    {
        return (int) Arrays.stream(depth).filter(d -> d > 0).count();
    }

    private static int[] calculateDepth(int contigLength, Collection<ViralReadAlignment> alignments)
    {
        int[] depth = new int[contigLength];
        for(ViralReadAlignment alignment : alignments)
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
