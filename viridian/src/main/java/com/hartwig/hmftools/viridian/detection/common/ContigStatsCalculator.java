package com.hartwig.hmftools.viridian.detection.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.stream.Collectors.toMap;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.align.AlignedInterval;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.reference.ViralContig;

// Statistics for every contig with an alignment, counting each read once by its best alignment there.
public class ContigStatsCalculator
{
    public static Map<ViralContig, ContigStats> calculate(
            Map<ViralContig, Map<ReadId, ViralReadAlignment>> alignmentsByContig,
            Map<ViralContig, Integer> originClippedReads)
    {
        // TODO: can this be a 1 statement stream operation?
        Map<ViralContig, Collection<ViralReadAlignment>> contigAlignments = new HashMap<>();
        alignmentsByContig.forEach((contig, readAlignments) -> contigAlignments.put(contig, readAlignments.values()));

        // A contig whose every alignment straddled the origin retains no read, but the drop is still worth reporting.
        originClippedReads.keySet().forEach(contig -> contigAlignments.computeIfAbsent(contig, k -> List.of()));

        return contigAlignments.entrySet().stream().collect(toMap(
                Map.Entry::getKey,
                entry -> createContigStats(
                        entry.getKey(), entry.getValue(), originClippedReads.getOrDefault(entry.getKey(), 0))));
    }

    private static ContigStats createContigStats(
            ViralContig contig, Collection<ViralReadAlignment> alignments, int originClippedReads)
    {
        int[] depth = calculateDepth(contig.length(), alignments);
        int coveredBases = (int) Arrays.stream(depth).filter(baseDepth -> baseDepth > 0).count();

        // Without retained reads there is no distribution to summarise, unlike depth which is zero everywhere.
        SummaryStats alignerScore = alignments.isEmpty()
                ? null : SummaryStats.from(alignments.stream().mapToInt(ViralReadAlignment::alignerScore).toArray());

        return new ContigStats(
                contig, alignments.size(), originClippedReads, coveredBases, SummaryStats.from(depth), alignerScore);
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
