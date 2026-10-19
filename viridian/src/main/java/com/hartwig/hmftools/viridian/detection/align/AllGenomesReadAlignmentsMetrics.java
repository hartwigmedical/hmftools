package com.hartwig.hmftools.viridian.detection.align;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.jetbrains.annotations.Nullable;

// Supporting information for the alignment of reads to all virus genomes.
public record AllGenomesReadAlignmentsMetrics(
        // Null only when no read aligned.
        @Nullable Double meanReadLength,
        // Reads with an alignment dropped for clipping over their virus genome's start or end.
        Map<VirusGenome, Integer> originClippedReads,
        // Per virus genome and per read, how many alignments the read had on the virus genome.
        Map<VirusGenome, List<Integer>> alignmentCountsByGenome
)
{
    public AllGenomesReadAlignmentsMetrics
    {
        boolean hasAlignments = !originClippedReads.isEmpty() || !alignmentCountsByGenome.isEmpty();
        if(hasAlignments == (meanReadLength == null))
        {
            throw new IllegalArgumentException(
                    "meanReadLength must be present when there are alignments: " + meanReadLength);
        }
        if(meanReadLength != null && meanReadLength <= 0)
        {
            throw new IllegalArgumentException("Invalid meanReadLength: " + meanReadLength);
        }
    }
}
