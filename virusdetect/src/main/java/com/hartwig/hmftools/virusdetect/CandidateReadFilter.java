package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateUnmapped;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

// Decides whether a read is a candidate for viral realignment.
public class CandidateReadFilter
{
    private final int mMinSoftClipBases;
    private final Set<String> mViralDecoyContigs;

    public CandidateReadFilter(int minSoftClipBases, Set<String> viralDecoyContigs)
    {
        mMinSoftClipBases = minSoftClipBases;
        mViralDecoyContigs = viralDecoyContigs;
    }

    public boolean isCandidate(SAMRecord record)
    {
        // Genuinely unaligned, so possibly viral. A read redux itself unmapped (UM tag) is instead host sequence
        // from a bad region, not genuinely unaligned, so it is excluded.
        if(record.getReadUnmappedFlag())
        {
            return !record.hasAttribute(UNMAP_ATTRIBUTE);
        }
        // The unmapped mate may be viral; this read anchors it.
        if(mateUnmapped(record))
        {
            return true;
        }
        // Mapped to a host-reference viral decoy contig.
        if(isViralDecoyContig(record.getReferenceName()))
        {
            return true;
        }
        // A long soft clip whose clipped bases were not placed elsewhere in the host may mark a viral junction.
        boolean hasSignificantClip = (leftSoftClipLength(record) >= mMinSoftClipBases || rightSoftClipLength(record) >= mMinSoftClipBases);
        if(hasSignificantClip && clippedBasesAreCandidate(record))
        {
            return true;
        }
        return false;
    }

    // A supplementary alignment places the clipped bases elsewhere in the host; if any lands on a non-viral contig the
    // clip is host sequence, not a viral junction. No supplementary means the clip may mark one.
    private boolean clippedBasesAreCandidate(SAMRecord record)
    {
        List<SupplementaryReadData> supplementaries = SupplementaryReadData.extractAlignments(record);
        if(supplementaries == null || supplementaries.isEmpty())
        {
            return true;
        }
        else
        {
            return supplementaries.stream().allMatch(supplementary -> isViralDecoyContig(supplementary.Chromosome));
        }
    }

    private boolean isViralDecoyContig(String contig)
    {
        return mViralDecoyContigs.contains(contig);
    }
}
