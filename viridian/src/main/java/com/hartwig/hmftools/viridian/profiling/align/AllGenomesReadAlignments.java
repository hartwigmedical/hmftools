package com.hartwig.hmftools.viridian.profiling.align;

import static java.util.stream.Collectors.toMap;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusReference;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Alignments of reads to (possibly multiple) virus genomes (via BWA-MEM -a mode).
// Alignments straddling a virus genome's origin are excluded here to avoid circular genome linearization artifacts.
public record AllGenomesReadAlignments(
        ViralReadAlignmentStore store,
        AllGenomesReadAlignmentsMetrics metrics
)
{
    public static AllGenomesReadAlignments from(List<ViralReadAlignment> alignments, @Nullable Double meanReadLength)
    {
        Map<VirusGenome, Set<ReadId>> originClippedReadsByGenome = new HashMap<>();
        Map<ReadGenome, List<ViralReadAlignment>> alignmentsByReadGenome = new LinkedHashMap<>();

        for(ViralReadAlignment alignment : alignments)
        {
            if(alignment.clipsOverGenomeEnd())
            {
                originClippedReadsByGenome
                        .computeIfAbsent(alignment.genome(), k -> new HashSet<>())
                        .add(alignment.readId());
            }
            else
            {
                alignmentsByReadGenome
                        .computeIfAbsent(new ReadGenome(alignment.readId(), alignment.genome()), k -> new ArrayList<>())
                        .add(alignment);
            }
        }

        List<ViralReadAlignment> retained = new ArrayList<>();
        Map<VirusGenome, List<Integer>> alignmentCountsByGenome = new HashMap<>();

        alignmentsByReadGenome.forEach((readGenome, readGenomeAlignments) ->
        {
            retained.add(readGenomeAlignments.stream().min(ViralReadAlignment.BEST_FIT_FIRST).orElseThrow());
            alignmentCountsByGenome
                    .computeIfAbsent(readGenome.genome(), k -> new ArrayList<>())
                    .add(readGenomeAlignments.size());
        });

        Map<VirusGenome, Integer> originClippedReads = originClippedReadsByGenome.entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> entry.getValue().size()));

        AllGenomesReadAlignmentsMetrics metrics =
                new AllGenomesReadAlignmentsMetrics(meanReadLength, originClippedReads, alignmentCountsByGenome);

        return new AllGenomesReadAlignments(new ViralReadAlignmentStore(retained), metrics);
    }

    public static AllGenomesReadAlignments load(String bamFile, VirusReference reference)
    {
        List<ViralReadAlignment> alignments = new ArrayList<>();
        long readLengthSum = 0;
        int readLengthCount = 0;

        try(SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile)))
        {
            for(SAMRecord record : reader)
            {
                if(record.getReadUnmappedFlag())
                {
                    continue;
                }
                alignments.add(ViralReadAlignment.from(record, reference));

                // Only primaries are significant for measuring read length.
                if(!record.isSecondaryAlignment() && !record.getSupplementaryAlignmentFlag() && record.getReadLength() > 0)
                {
                    readLengthSum += record.getReadLength();
                    ++readLengthCount;
                }
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to read aligned BAM", e);
        }

        Double meanReadLength = readLengthCount > 0 ? (double) readLengthSum / readLengthCount : null;
        return from(alignments, meanReadLength);
    }

    private record ReadGenome(
            ReadId readId,
            VirusGenome genome
    )
    {
    }
}
