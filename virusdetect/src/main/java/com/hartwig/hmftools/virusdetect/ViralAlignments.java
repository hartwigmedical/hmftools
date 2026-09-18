package com.hartwig.hmftools.virusdetect;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// All viral alignments read from the viral reference alignment BAM.
// Each read may have multiple alignments (BWA-MEM -a mode).
// Alignments straddling a contig's origin are excluded here to avoid linearization artifacts.
public record ViralAlignments(
        List<ViralAlignment> alignments,
        double meanReadLength,
        // Alignments excluded for clipping over their contig's start or end.
        Map<ViralContig, Integer> originClippedReads,
        // Distinct reads with at least one alignment to any contig of the oncology group.
        Map<String, Integer> readCountsByOncologyGroup
)
{
    public ViralAlignments
    {
        if(alignments.stream().anyMatch(ViralAlignment::clipsOverContigEnd))
        {
            throw new IllegalArgumentException("Origin-straddling alignments must be excluded");
        }
    }

    public static ViralAlignments from(List<ViralAlignment> alignments, double meanReadLength)
    {
        List<ViralAlignment> retained = new ArrayList<>();
        Map<ViralContig, Integer> originClippedReads = new HashMap<>();
        Map<String, Set<String>> readNamesByOncologyGroup = new HashMap<>();

        for(ViralAlignment alignment : alignments)
        {
            if(alignment.clipsOverContigEnd())
            {
                originClippedReads.merge(alignment.contig(), 1, Integer::sum);
                continue;
            }

            retained.add(alignment);
            readNamesByOncologyGroup
                    .computeIfAbsent(alignment.contig().oncologyGroup(), oncologyGroup -> new HashSet<>())
                    .add(alignment.readName());
        }

        Map<String, Integer> readCounts = new HashMap<>();
        readNamesByOncologyGroup.forEach((oncologyGroup, readNames) -> readCounts.put(oncologyGroup, readNames.size()));

        return new ViralAlignments(retained, meanReadLength, originClippedReads, readCounts);
    }

    public static ViralAlignments load(String bamFile, ViralReference reference)
    {
        List<ViralAlignment> alignments = new ArrayList<>();
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
                alignments.add(ViralAlignment.from(record, reference));

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
            throw new RuntimeException("failed to read aligned BAM", e);
        }

        double meanReadLength = readLengthCount > 0 ? (double) readLengthSum / readLengthCount : 0.0;
        return from(alignments, meanReadLength);
    }
}
