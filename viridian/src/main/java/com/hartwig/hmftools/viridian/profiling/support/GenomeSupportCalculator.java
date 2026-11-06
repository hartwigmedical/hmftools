package com.hartwig.hmftools.viridian.profiling.support;

import static java.util.stream.Collectors.toMap;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.profiling.align.AllGenomesReadAlignments;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStatsCalculator;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Assembles per-genome read support information, particularly read attribution counts used for representative selection.
public class GenomeSupportCalculator
{
    private final double mCorrectBaseProbability;

    public GenomeSupportCalculator(double correctBaseProbability)
    {
        mCorrectBaseProbability = correctBaseProbability;
    }

    public List<GenomeSupport> calculate(AllGenomesReadAlignments alignments)
    {
        Map<VirusGenome, GenomeStats> genomeStats = GenomeStatsCalculator.calculate(
                alignments.genomes(), alignments.store().byGenome(), alignments.originClippedReads());
        Collection<Map<VirusGenome, ViralReadAlignment>> reads = alignments.store().byRead().values();
        Map<VirusGenome, Double> readVotes = calculateReadVotes(reads);
        Map<VirusGenome, Integer> topMatchReads = calculateTopMatchReads(reads);
        Map<VirusGenome, List<Integer>> alignmentCounts = alignments.alignmentCountsByGenome();

        return genomeStats.values().stream()
                .map(stats -> createGenomeSupport(
                        stats, readVotes.getOrDefault(stats.genome(), 0.0),
                        topMatchReads.getOrDefault(stats.genome(), 0),
                        alignmentCounts.getOrDefault(stats.genome(), List.of())))
                .toList();
    }

    // "Read votes":
    // Each read gives total 1 vote, split between its aligned genomes according to edit distance.
    // As edit distance to the genome increases, the vote to that genome decreases (distributed to the other aligned genomes).
    // Example:
    //   If a read aligned to 2 genomes, one perfectly and one with a mismatch, then the votes would be 0.8 and 0.2.
    // Downside is that with many similar genomes, the vote counts are diluted, so all genomes appear to be equally good.
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

    // "Top match reads":
    // For each genome, counts the reads whose overall best alignment is to that genome.
    // A read with tied best alignments counts for each tied genome, so unlike read votes, the count does not dilute
    // with many similar genomes.
    private static Map<VirusGenome, Integer> calculateTopMatchReads(
            Collection<Map<VirusGenome, ViralReadAlignment>> reads)
    {
        Map<VirusGenome, Integer> topMatchReads = new HashMap<>();
        for(Map<VirusGenome, ViralReadAlignment> readAlignments : reads)
        {
            // Find the lowest divergence, i.e. best alignment.
            int minDivergence = readAlignments.values().stream()
                    .mapToInt(ViralReadAlignment::divergence).min().orElseThrow();
            // Get all the genomes which align equally well as the best.
            readAlignments.entrySet().stream()
                    .filter(entry -> entry.getValue().divergence() == minDivergence)
                    .forEach(entry -> topMatchReads.merge(entry.getKey(), 1, Integer::sum));
        }
        return topMatchReads;
    }

    private static GenomeSupport createGenomeSupport(
            GenomeStats stats, double readVotes, int topMatchReads, List<Integer> alignmentCounts)
    {
        int multiAlignReads = (int) alignmentCounts.stream().filter(count -> count > 1).count();

        // Without retained reads there is no distribution to summarise.
        SummaryStats alignmentsPerRead = alignmentCounts.isEmpty() ? null : SummaryStats.from(alignmentCounts);

        return new GenomeSupport(stats, multiAlignReads, alignmentsPerRead, readVotes, topMatchReads);
    }
}
