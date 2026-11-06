package com.hartwig.hmftools.viridian.profiling.common;

import static java.util.function.Function.identity;
import static java.util.stream.Collectors.toMap;

import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.viridian.profiling.align.MatchedInterval;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Statistics for every virus genome with an alignment, counting each read once by its best alignment there.
public class GenomeStatsCalculator
{
    public static Map<VirusGenome, GenomeStats> calculate(
            Set<VirusGenome> genomes, Map<VirusGenome, Map<ReadId, ViralReadAlignment>> alignmentsByGenome,
            Map<VirusGenome, Integer> originClippedReads)
    {
        return genomes.stream()
                .collect(toMap(
                        identity(),
                        genome -> createGenomeStats(
                                genome,
                                alignmentsByGenome.getOrDefault(genome, Map.of()).values(),
                                originClippedReads.getOrDefault(genome, 0))));
    }

    private static GenomeStats createGenomeStats(
            VirusGenome genome, Collection<ViralReadAlignment> alignments, int originClippedReads)
    {
        int[] depth = calculateDepth(genome.length(), alignments);
        int coveredBases = (int) Arrays.stream(depth).filter(baseDepth -> baseDepth > 0).count();

        // Without retained reads there is no distribution to summarise, unlike depth which is zero everywhere.
        SummaryStats alignerScore = alignments.isEmpty()
                ? null : SummaryStats.from(alignments.stream().mapToInt(ViralReadAlignment::alignerScore).toArray());

        return new GenomeStats(
                genome, alignments.size(), originClippedReads, coveredBases, SummaryStats.from(depth), alignerScore);
    }

    private static int[] calculateDepth(int genomeLength, Collection<ViralReadAlignment> alignments)
    {
        int[] depth = new int[genomeLength];
        for(ViralReadAlignment alignment : alignments)
        {
            for(MatchedInterval block : alignment.matchedIntervals())
            {
                int start = block.referenceStart() - 1;   // intervals are 1-based
                int end = start + block.length();
                for(int position = start; position < end; ++position)
                {
                    ++depth[position];
                }
            }
        }
        return depth;
    }
}
