package com.hartwig.hmftools.virusdetect;

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;

import java.util.List;
import java.util.Map;

// One read's viral alignments, reduced to its best alignment on each contig it hits.
public record ReadAlignments(
        String readName,
        Map<ViralContig, ReadContigAlignment> hits
)
{
    public static ReadAlignments from(String readName, List<ViralAlignment> alignments)
    {
        Map<ViralContig, ReadContigAlignment> hits = alignments.stream()
                .collect(groupingBy(ViralAlignment::contig))
                .entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> ReadContigAlignment.from(entry.getValue())));

        return new ReadAlignments(readName, hits);
    }

    // The fewest divergent bases across the contigs the read hits. Offsets vote weighting for numeric stability.
    public int minDivergence()
    {
        return hits.values().stream().mapToInt(ReadContigAlignment::divergence).min().orElseThrow();
    }
}
