package com.hartwig.hmftools.viridian.detection.support;

import static java.util.Objects.requireNonNull;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.align.AllContigsReadAlignments;
import com.hartwig.hmftools.viridian.detection.align.AllContigsReadAlignmentsMetrics;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.detection.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Takes per-contig statistics and determines candidate virus genomes.
// Computes:
//   - Baseline support status (is virus maybe here or is it noise?);
//   - Read voting used for representative contig selection.
public class ContigSupportCalculator
{
    private final double mCorrectBaseProbability;

    public ContigSupportCalculator(double correctBaseProbability)
    {
        mCorrectBaseProbability = correctBaseProbability;
    }

    public List<ContigSupport> calculate(AllContigsReadAlignments alignments, Map<VirusGenome, ContigStats> contigStats)
    {
        Map<VirusGenome, Double> readVotes = calculateReadVotes(alignments.store().byRead().values());
        AllContigsReadAlignmentsMetrics metrics = alignments.metrics();

        return contigStats.values().stream()
                .collect(groupingBy(stats -> stats.contig().oncologyGroup()))
                .values().stream()
                .flatMap(group -> createGroupContigSupports(
                        group, readVotes, metrics.alignmentCountsByContig(),
                        // Contig stats exist only for contigs with alignments, so the mean read length is present.
                        requireNonNull(metrics.meanReadLength())).stream())
                .toList();
    }

    // A read's vote splits across the contigs it hits by how well each explains it: every extra base a contig fails to
    // explain multiplies its share by the (pessimistic) chance a base is right, so the closest contig wins most.
    private Map<VirusGenome, Double> calculateReadVotes(Collection<Map<VirusGenome, ViralReadAlignment>> reads)
    {
        Map<VirusGenome, Double> votes = new HashMap<>();
        for(Map<VirusGenome, ViralReadAlignment> readAlignments : reads)
        {
            // Offsets vote weighting for numeric stability.
            int minDivergence = readAlignments.values().stream()
                    .mapToInt(ViralReadAlignment::divergence).min().orElseThrow();
            Map<VirusGenome, Double> contigWeights = readAlignments.entrySet().stream()
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
            List<ContigStats> group, Map<VirusGenome, Double> readVotes, Map<VirusGenome, List<Integer>> alignmentCounts,
            double meanReadLength)
    {
        Map<VirusGenome, ContigStats> groupStats = group.stream().collect(toMap(ContigStats::contig, stats -> stats));
        Map<VirusGenome, ContigFilterStatus> filterStatuses = ContigSupportFilter.statuses(groupStats, readVotes, meanReadLength);

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
