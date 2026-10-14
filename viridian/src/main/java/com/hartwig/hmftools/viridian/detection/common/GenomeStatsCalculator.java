package com.hartwig.hmftools.viridian.detection.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.function.Function.identity;
import static java.util.stream.Collectors.toMap;

import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.stream.Stream;

import com.hartwig.hmftools.viridian.detection.align.AlignedInterval;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Statistics for every virus genome with an alignment, counting each read once by its best alignment there.
public class GenomeStatsCalculator
{
    public static Map<VirusGenome, GenomeStats> calculate(
            Map<VirusGenome, Map<ReadId, ViralReadAlignment>> alignmentsByGenome,
            Map<VirusGenome, Integer> originClippedReads)
    {
        // A virus genome whose every alignment straddled the origin retains no read, but the drop is still worth
        // reporting.
        return Stream.concat(alignmentsByGenome.keySet().stream(), originClippedReads.keySet().stream())
                .distinct()
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
            for(AlignedInterval block : alignment.alignedIntervals())
            {
                int start = max(0, block.referenceStart() - 1);   // intervals are 1-based
                int end = min(genomeLength, block.referenceStart() - 1 + block.length());
                for(int position = start; position < end; ++position)
                {
                    ++depth[position];
                }
            }
        }
        return depth;
    }
}
