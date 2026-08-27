package com.hartwig.hmftools.virusdetect;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateUnmapped;

import java.util.Set;

import htsjdk.samtools.SAMRecord;

// Decides whether a read is a candidate for viral realignment: fully unmapped, one mate unmapped,
// long soft clip, or mapped to a host-reference decoy contig.
public class CandidateReadFilter
{
    private final int mMinSoftClipBases;
    private final Set<String> mDecoyContigs;

    public CandidateReadFilter(int minSoftClipBases, Set<String> decoyContigs)
    {
        mMinSoftClipBases = minSoftClipBases;
        mDecoyContigs = decoyContigs;
    }

    public boolean isCandidate(SAMRecord record)
    {
        if(record.getReadUnmappedFlag())
        {
            return true;
        }
        if(mateUnmapped(record))
        {
            return true;
        }
        if(max(leftSoftClipLength(record), rightSoftClipLength(record)) >= mMinSoftClipBases)
        {
            return true;
        }
        return mDecoyContigs.contains(record.getReferenceName());
    }
}
