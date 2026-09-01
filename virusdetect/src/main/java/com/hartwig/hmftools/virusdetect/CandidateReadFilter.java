package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateUnmapped;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

// Decides whether a read is a candidate for viral realignment: fully unmapped, one mate unmapped, mapped to a
// host-reference decoy contig, or long soft clip whose clipped bases are not already placed elsewhere in the host.
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
        if((leftSoftClipLength(record) >= mMinSoftClipBases || rightSoftClipLength(record) >= mMinSoftClipBases)
                && clippedBasesNotInHost(record))
        {
            return true;
        }
        return mDecoyContigs.contains(record.getReferenceName());
    }

    // A supplementary alignment means the aligner placed the clipped bases elsewhere in the host reference; if any
    // lands on a non-viral contig the clip is explained by host sequence, not a viral junction. With no supplementary
    // the clipped bases were not placed in the host, so the clip may mark a viral junction.
    private boolean clippedBasesNotInHost(SAMRecord record)
    {
        List<SupplementaryReadData> supplementaries = SupplementaryReadData.extractAlignments(record);
        if(supplementaries == null)
        {
            return true;
        }
        return supplementaries.stream().allMatch(supplementary -> mDecoyContigs.contains(supplementary.Chromosome));
    }
}
