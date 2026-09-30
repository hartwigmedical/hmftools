package com.hartwig.hmftools.viridian.detection.support;

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.READ_VOTE_CORRECT_BASE_PROBABILITY;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.align.AllAlignments;
import com.hartwig.hmftools.viridian.detection.align.AllAlignmentsMetrics;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.detection.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.ViralContig;

// Turns per-contig statistics into a verdict on the contig being present in the sample.
// Computes:
//   - Baseline presence status;
//   - Read voting used for representative contig selection.
public class ContigSupportCalculator
{
    private final double mCorrectBaseProbability;

    public ContigSupportCalculator(double correctBaseProbability)
    {
        mCorrectBaseProbability = correctBaseProbability;
    }

    public ContigSupportCalculator()
    {
        this(READ_VOTE_CORRECT_BASE_PROBABILITY);
    }

    // TODO rename "calculate"
    public List<ContigSupport> compute(AllAlignments allAlignments, Map<ViralContig, ContigStats> contigStats)
    {
        Map<ViralContig, Double> readVotes = calculateReadVotes(allAlignments.alignments().byRead().values());
        AllAlignmentsMetrics metrics = allAlignments.metrics();

        return contigStats.values().stream()
                .collect(groupingBy(stats -> stats.contig().oncologyGroup()))
                .values().stream()
                .flatMap(group -> createGroupContigSupports(
                        group, readVotes, metrics.alignmentCountsByContig(), metrics.meanReadLength()).stream())
                .toList();
    }

    // A read's vote splits across the contigs it hits by how well each explains it: every extra base a contig fails to
    // explain multiplies its share by the (pessimistic) chance a base is right, so the closest contig wins most.
    private Map<ViralContig, Double> calculateReadVotes(Collection<Map<ViralContig, ViralReadAlignment>> reads)
    {
        Map<ViralContig, Double> votes = new HashMap<>();
        for(Map<ViralContig, ViralReadAlignment> readAlignments : reads)
        {
            // Offsets vote weighting for numeric stability.
            int minDivergence = readAlignments.values().stream()
                    .mapToInt(ViralReadAlignment::divergence).min().orElseThrow();
            Map<ViralContig, Double> contigWeights = readAlignments.entrySet().stream()
                    .collect(toMap(Map.Entry::getKey, entry -> voteWeight(entry.getValue().divergence(), minDivergence)));

            double totalWeight = contigWeights.values().stream().mapToDouble(Double::doubleValue).sum();
            contigWeights.forEach((contig, weight) -> votes.merge(contig, weight / totalWeight, Double::sum));
        }
        return votes;
    }

    private double voteWeight(int divergence, int minDivergence)
    {
        return Math.pow(mCorrectBaseProbability, divergence - minDivergence);
    }

    private static List<ContigSupport> createGroupContigSupports(
            List<ContigStats> group, Map<ViralContig, Double> readVotes, Map<ViralContig, List<Integer>> alignmentCounts,
            double meanReadLength)
    {
        Map<ViralContig, ContigStats> groupStats = group.stream().collect(toMap(ContigStats::contig, stats -> stats));
        Map<ViralContig, ContigFilterStatus> filterStatuses = ContigSupportFilter.statuses(groupStats, readVotes, meanReadLength);

        return group.stream()
                .map(stats -> createContigSupport(
                        stats, filterStatuses.get(stats.contig()), readVotes.getOrDefault(stats.contig(), 0.0),
                        alignmentCounts.getOrDefault(stats.contig(), List.of())))
                .toList();
    }

    private static ContigSupport createContigSupport(
            ContigStats stats, ContigFilterStatus filterStatus, double readVotes, List<Integer> alignmentCounts)
    {
        int multiAlignReads = (int) alignmentCounts.stream().filter(count -> count > 1).count();

        // Without retained reads there is no distribution to summarise.
        SummaryStats alignmentsPerRead = alignmentCounts.isEmpty() ? null : SummaryStats.from(alignmentCounts);

        return new ContigSupport(stats, filterStatus, multiAlignReads, alignmentsPerRead, readVotes);
    }
}
