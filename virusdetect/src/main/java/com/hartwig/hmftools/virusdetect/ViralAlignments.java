package com.hartwig.hmftools.virusdetect;

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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// All viral alignments read from the viral reference alignment BAM.
// Each read may have multiple alignments (BWA-MEM -a mode).
// Alignments straddling a contig's origin are excluded here to avoid linearization artifacts.
public record ViralAlignments(
        List<ReadAlignments> reads,
        double meanReadLength,
        // Reads with an alignment dropped for clipping over their contig's start or end.
        Map<ViralContig, Integer> originClippedReads,
        // Distinct reads with at least one alignment to any contig of the oncology group.
        Map<OncologyGroup, Integer> readCountsByOncologyGroup
)
{
    public ViralAlignments
    {
        if(!reads.isEmpty() && meanReadLength <= 0)
        {
            throw new IllegalArgumentException("Invalid mean read length: " + meanReadLength);
        }
    }

    public static ViralAlignments from(List<ViralAlignment> alignments, double meanReadLength)
    {
        Map<ViralContig, Set<String>> originClippedReadsByContig = new HashMap<>();
        Map<OncologyGroup, Set<String>> readsByOncologyGroup = new HashMap<>();
        Map<String, List<ViralAlignment>> alignmentsByRead = new LinkedHashMap<>();

        for(ViralAlignment alignment : alignments)
        {
            if(alignment.clipsOverContigEnd())
            {
                Set<String> contigReads = originClippedReadsByContig.computeIfAbsent(alignment.contig(), k -> new HashSet<>());
                contigReads.add(alignment.readName());
            }
            else{
                List<ViralAlignment> readAlignments = alignmentsByRead.computeIfAbsent(alignment.readName(), k -> new ArrayList<>());
                readAlignments.add(alignment);

                Set<String> oncologyGroupReads = readsByOncologyGroup.computeIfAbsent(alignment.contig().oncologyGroup(), k -> new HashSet<>());
                oncologyGroupReads.add(alignment.readName());
            }
        }

        List<ReadAlignments> reads = alignmentsByRead.entrySet().stream()
                .map(entry -> ReadAlignments.from(entry.getKey(), entry.getValue()))
                .toList();

        Map<ViralContig, Integer> originClippedReadCounts = originClippedReadsByContig.entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> entry.getValue().size()));

        Map<OncologyGroup, Integer> oncologyGroupReadCounts = readsByOncologyGroup.entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> entry.getValue().size()));

        return new ViralAlignments(reads, meanReadLength, originClippedReadCounts, oncologyGroupReadCounts);
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
            throw new RuntimeException("Failed to read aligned BAM", e);
        }

        double meanReadLength = readLengthCount > 0 ? (double) readLengthSum / readLengthCount : 0.0;
        return from(alignments, meanReadLength);
    }
}
