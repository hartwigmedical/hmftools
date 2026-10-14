package com.hartwig.hmftools.viridian.detection.support;

import static java.util.Objects.requireNonNull;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.align.AllGenomesReadAlignments;
import com.hartwig.hmftools.viridian.detection.align.AllGenomesReadAlignmentsMetrics;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.detection.common.GenomeStats;
import com.hartwig.hmftools.viridian.detection.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Takes per-virus-genome statistics and determines candidate virus genomes.
// Computes:
//   - Baseline support status (is virus maybe here or is it noise?);
//   - Read voting used for representative selection.
public class GenomeSupportCalculator
{
    private final double mCorrectBaseProbability;

    public GenomeSupportCalculator(double correctBaseProbability)
    {
        mCorrectBaseProbability = correctBaseProbability;
    }

    public List<GenomeSupport> calculate(AllGenomesReadAlignments alignments, Map<VirusGenome, GenomeStats> genomeStats)
    {
        Map<VirusGenome, Double> readVotes = calculateReadVotes(alignments.store().byRead().values());
        AllGenomesReadAlignmentsMetrics metrics = alignments.metrics();

        return genomeStats.values().stream()
                .collect(groupingBy(stats -> stats.genome().oncologyGroup()))
                .values().stream()
                .flatMap(group -> createGroupGenomeSupports(
                        group, readVotes, metrics.alignmentCountsByGenome(),
                        // Stats exist only for virus genomes with alignments, so the mean read length is present.
                        requireNonNull(metrics.meanReadLength())).stream())
                .toList();
    }

    // A read's vote splits across the virus genomes it hits by how well each explains it: every extra base a genome
    // fails to explain multiplies its share by the (pessimistic) chance a base is right, so the closest genome wins
    // most.
    private Map<VirusGenome, Double> calculateReadVotes(Collection<Map<VirusGenome, ViralReadAlignment>> reads)
    {
        Map<VirusGenome, Double> votes = new HashMap<>();
        for(Map<VirusGenome, ViralReadAlignment> readAlignments : reads)
        {
            // Offsets vote weighting for numeric stability.
            int minDivergence = readAlignments.values().stream()
                    .mapToInt(ViralReadAlignment::divergence).min().orElseThrow();
            Map<VirusGenome, Double> genomeWeights = readAlignments.entrySet().stream()
                    .collect(toMap(Map.Entry::getKey, entry -> voteWeight(entry.getValue().divergence(), minDivergence)));

            double totalWeight = genomeWeights.values().stream().mapToDouble(Double::doubleValue).sum();
            genomeWeights.forEach((genome, weight) -> votes.merge(genome, weight / totalWeight, Double::sum));
        }
        return votes;
    }

    private double voteWeight(int divergence, int minDivergence)
    {
        return Math.pow(mCorrectBaseProbability, divergence - minDivergence);
    }

    private static List<GenomeSupport> createGroupGenomeSupports(
            List<GenomeStats> group, Map<VirusGenome, Double> readVotes, Map<VirusGenome, List<Integer>> alignmentCounts,
            double meanReadLength)
    {
        Map<VirusGenome, GenomeStats> groupStats = group.stream().collect(toMap(GenomeStats::genome, stats -> stats));
        Map<VirusGenome, GenomeFilterStatus> filterStatuses = GenomeSupportFilter.statuses(groupStats, readVotes, meanReadLength);

        return group.stream()
                .map(stats -> createGenomeSupport(
                        stats, filterStatuses.get(stats.genome()), readVotes.getOrDefault(stats.genome(), 0.0),
                        alignmentCounts.getOrDefault(stats.genome(), List.of())))
                .toList();
    }

    private static GenomeSupport createGenomeSupport(
            GenomeStats stats, GenomeFilterStatus filterStatus, double readVotes, List<Integer> alignmentCounts)
    {
        int multiAlignReads = (int) alignmentCounts.stream().filter(count -> count > 1).count();

        // Without retained reads there is no distribution to summarise.
        SummaryStats alignmentsPerRead = alignmentCounts.isEmpty() ? null : SummaryStats.from(alignmentCounts);

        return new GenomeSupport(stats, filterStatus, multiAlignReads, alignmentsPerRead, readVotes);
    }
}
