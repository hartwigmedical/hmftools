package com.hartwig.hmftools.virusdetect;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// All viral alignments read from the viral reference alignment BAM.
// Each read may have multiple alignments (BWA-MEM -a mode)
public record ViralAlignments(
        List<ViralAlignment> alignments,
        double meanReadLength
)
{
    public static ViralAlignments load(String bamFile)
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
                alignments.add(ViralAlignment.from(record));

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
        return new ViralAlignments(alignments, meanReadLength);
    }
}
