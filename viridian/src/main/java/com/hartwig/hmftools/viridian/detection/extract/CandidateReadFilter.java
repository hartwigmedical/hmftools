package com.hartwig.hmftools.viridian.detection.extract;

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
    private final Set<String> mRefVirusContigs;

    // TODO: add 0-arg construct that uses constants

    public CandidateReadFilter(int minSoftClipBases, Set<String> refVirusContigs)
    {
        mMinSoftClipBases = minSoftClipBases;
        mRefVirusContigs = refVirusContigs;
    }

    public boolean isCandidate(SAMRecord record)
    {
        // Only care about primaries because only they have the read sequence.
        // BamSlicer filters these anyway, but check here just in case.
        if(record.isSecondaryOrSupplementary())
        {
            return false;
        }
        // REDUX duplicate fragment. Doesn't provide additional support.
        // BamSlicer filters these anyway, but check here just in case.
        if(record.getDuplicateReadFlag())
        {
            return false;
        }
        // Mapped to a virus decoy contig, or an unmapped read placed on one by its mapped mate. Either way the fragment
        // touches a virus, so keep it. Checked before the redux-unmapped exclusion, since a redux-unmapped read sitting
        // on a decoy is still viral evidence.
        if(isVirusDecoyContig(record.getReferenceName()))
        {
            return true;
        }
        // The mate maps to a virus decoy contig, so this read anchors a host<->virus fragment (e.g. an integration junction).
        if(mateMappedToVirusDecoy(record))
        {
            return true;
        }
        // Genuinely unaligned, so possibly viral. A read redux itself unmapped (UM tag) is instead host sequence
        // from a bad region, not genuinely unaligned, so it is not a candidate on its own account.
        if(record.getReadUnmappedFlag())
        {
            return !record.hasAttribute(UNMAP_ATTRIBUTE);
        }
        // The unmapped mate may be viral, so take this read for integration support.
        if(mateUnmapped(record))
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

    private boolean clippedBasesAreCandidate(SAMRecord record)
    {
        List<SupplementaryReadData> supplementaries = SupplementaryReadData.extractAlignments(record);
        if(supplementaries == null || supplementaries.isEmpty())
        {
            return true;
        }
        else
        {
            // We assume that if ANY supplementary is possibly host genome, then it's not viral, even if another supplementary maps to a
            // virus contig.
            return supplementaries.stream().allMatch(supplementary -> isVirusDecoyContig(supplementary.Chromosome));
        }
    }

    private boolean mateMappedToVirusDecoy(SAMRecord record)
    {
        return record.getReadPairedFlag() && !record.getMateUnmappedFlag() && isVirusDecoyContig(record.getMateReferenceName());
    }

    private boolean isVirusDecoyContig(String contig)
    {
        return mRefVirusContigs.contains(contig);
    }
}
