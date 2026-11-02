package com.hartwig.hmftools.viridian.profiling.support;

import static java.util.stream.Collectors.toMap;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.profiling.align.AllGenomesReadAlignments;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Assembles per-genome read support information and computes read voting used for representative selection.
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
        Map<VirusGenome, List<Integer>> alignmentCounts = alignments.metrics().alignmentCountsByGenome();

        return genomeStats.values().stream()
                .map(stats -> createGenomeSupport(
                        stats, readVotes.getOrDefault(stats.genome(), 0.0),
                        alignmentCounts.getOrDefault(stats.genome(), List.of())))
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

    private static GenomeSupport createGenomeSupport(GenomeStats stats, double readVotes, List<Integer> alignmentCounts)
    {
        int multiAlignReads = (int) alignmentCounts.stream().filter(count -> count > 1).count();

        // Without retained reads there is no distribution to summarise.
        SummaryStats alignmentsPerRead = alignmentCounts.isEmpty() ? null : SummaryStats.from(alignmentCounts);

        return new GenomeSupport(stats, multiAlignReads, alignmentsPerRead, readVotes);
    }
}
