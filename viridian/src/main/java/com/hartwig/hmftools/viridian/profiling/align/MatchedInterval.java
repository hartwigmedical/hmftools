package com.hartwig.hmftools.viridian.profiling.align;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE;
import static com.hartwig.hmftools.viridian.common.Utils.getRequiredSamTag;
import static com.hartwig.hmftools.viridian.common.Utils.parseMdTagMismatchPositions;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;

// A contiguous run of reference bases matched in an alignment.
public record MatchedInterval(
        // 1-based inclusive.
        int referenceStart,
        int length
)
{
    public MatchedInterval
    {
        if(referenceStart < 1)
        {
            throw new IllegalArgumentException("Invalid reference start: " + referenceStart);
        }
        if(length < 1)
        {
            throw new IllegalArgumentException("Invalid interval length: " + length);
        }
    }

    public static List<MatchedInterval> fromSamRecord(SAMRecord record)
    {
        String mdTag = getRequiredSamTag(record, MISMATCHES_AND_DELETIONS_ATTRIBUTE, SAMRecord::getStringAttribute);
        List<Integer> mismatchPositions = parseMdTagMismatchPositions(record.getAlignmentStart(), mdTag);

        List<MatchedInterval> intervals = new ArrayList<>();
        int nextMismatch = 0;
        for(AlignmentBlock block : record.getAlignmentBlocks())
        {
            int start = block.getReferenceStart();
            int end = start + block.getLength();
            for(; nextMismatch < mismatchPositions.size() && mismatchPositions.get(nextMismatch) < end; ++nextMismatch)
            {
                int mismatch = mismatchPositions.get(nextMismatch);
                if(mismatch < start)
                {
                    throw inconsistentMdTag(record);
                }
                if(mismatch > start)
                {
                    intervals.add(new MatchedInterval(start, mismatch - start));
                }
                start = mismatch + 1;
            }
            if(end > start)
            {
                intervals.add(new MatchedInterval(start, end - start));
            }
        }
        if(nextMismatch < mismatchPositions.size())
        {
            throw inconsistentMdTag(record);
        }
        return intervals;
    }

    private static IllegalStateException inconsistentMdTag(SAMRecord record)
    {
        return new IllegalStateException(
                String.format("aligned read MD tag disagrees with its CIGAR on contig: %s", record.getReferenceName()));
    }
}
