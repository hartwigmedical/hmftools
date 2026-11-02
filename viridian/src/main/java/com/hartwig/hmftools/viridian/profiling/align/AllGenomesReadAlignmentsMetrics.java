package com.hartwig.hmftools.viridian.profiling.align;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Supporting information for the alignment of reads to all virus genomes.
public record AllGenomesReadAlignmentsMetrics(
        // Reads with an alignment dropped for clipping over their virus genome's start or end.
        Map<VirusGenome, Integer> originClippedReads,
        // Per virus genome and per read, how many alignments the read had on the virus genome.
        Map<VirusGenome, List<Integer>> alignmentCountsByGenome
)
{
}
