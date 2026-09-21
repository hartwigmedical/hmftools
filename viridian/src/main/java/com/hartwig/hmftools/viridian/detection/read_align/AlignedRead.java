package com.hartwig.hmftools.viridian.detection.read_align;

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.reference.ViralContig;

// One read's viral alignments, reduced to its best alignment on each contig with alignments.
public record AlignedRead(
        String readName,
        Map<ViralContig, ReadContigAlignment> hits
)
{
    public static AlignedRead from(String readName, List<ViralReadAlignment> alignments)
    {
        Map<ViralContig, ReadContigAlignment> hits = alignments.stream()
                .collect(groupingBy(ViralReadAlignment::contig))
                .entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> ReadContigAlignment.from(entry.getValue())));

        return new AlignedRead(readName, hits);
    }

    // The fewest divergent bases across the contigs the read hits. Offsets vote weighting for numeric stability.
    public int minDivergence()
    {
        return hits.values().stream().mapToInt(ReadContigAlignment::divergence).min().orElseThrow();
    }
}
