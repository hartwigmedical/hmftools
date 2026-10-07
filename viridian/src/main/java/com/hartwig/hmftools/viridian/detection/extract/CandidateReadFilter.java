package com.hartwig.hmftools.viridian.detection.extract;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAPP_COORDS_DELIM;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

// Decides whether a read is a candidate for viral realignment.
public class CandidateReadFilter
{
    private final CandidateReadSource mSource;
    private final int mMinSoftClipBases;
    private final Set<String> mRefVirusContigs;
    @Nullable
    private final ViralKmerIndex mKmerIndex;

    public CandidateReadFilter(
            CandidateReadSource source, int minSoftClipBases, Set<String> refVirusContigs,
            @Nullable ViralKmerIndex kmerIndex)
    {
        mSource = source;
        mMinSoftClipBases = minSoftClipBases;
        mRefVirusContigs = refVirusContigs;
        mKmerIndex = kmerIndex;
    }

    public boolean isCandidateContig(String contig)
    {
        return switch(mSource)
        {
            case ALL -> true;
            case FULLY_UNMAPPED_AND_VIRUS -> isVirusDecoyContig(contig);
        };
    }

    public boolean isCandidateRecord(SAMRecord record)
    {
        return isStructuralCandidate(record) && hasViralKmer(record);
    }

    private boolean isStructuralCandidate(SAMRecord record)
    {
        if(!isCandidateRecordType(record))
        {
            return false;
        }

        boolean isReduxUnmappedFromHost = false;
        // Consider evidence that the read is viral from the read itself.
        if(record.getReadUnmappedFlag())
        {
            if(isReduxUnmapped(record))
            {
                if(isReduxUnmappedFromVirusDecoy(record))
                {
                    // REDUX unmapped but the original contig was virus decoy: definitely viral fragment.
                    return true;
                }
                else
                {
                    // REDUX unmapped and the original contig was regular human genome: no direct viral evidence for the read; fall through
                    // to mate consideration.
                    isReduxUnmappedFromHost = true;
                }
            }
            else
            {
                // Genuinely unaligned read: maybe viral fragment.
                return true;
            }
        }
        else
        {
            if(isVirusDecoyContig(record.getReferenceName()))
            {
                // Mapped to a virus decoy contig in the ref genome: definitely viral fragment.
                return true;
            }
            else if(hasCandidateSoftClip(record))
            {
                // Mapped but clipped so could be a viral integration: maybe viral fragment.
                return true;
            }
            else
            {
                // Mapped to regular human genome: no direct viral evidence for the read; fall through to mate consideration.
            }
        }

        // Now consider evidence that the read is viral from the mate.
        if(record.getReadPairedFlag())
        {
            if(record.getMateUnmappedFlag())
            {
                if(isReduxUnmappedFromHost)
                {
                    // Annoying case: When the read was REDUX unmapped and its mate is unmapped too, the mate was almost always REDUX
                    // unmapped as well, so it isn't viral evidence for the read.
                    // A genuinely unaligned mate is still accepted when the mate record itself is checked.
                    return false;
                }
                else
                {
                    // Mate genuinely unaligned: maybe viral fragment.
                    // OR the mate is unmapped by REDUX, but there is no way to know, so have to be conservative.
                    return true;
                }
            }
            else
            {
                if(isVirusDecoyContig(record.getMateReferenceName()))
                {
                    // Mate mapped to a virus decoy contig in the ref genome: definitely viral fragment.
                    return true;
                }
                else
                {
                    // TODO: can read mate cigar attribute to look at soft clip of mate?

                    // Mapped to regular human genome: no viral evidence for the mate.
                }
            }
        }

        // No evidence remaining: probably not viral fragment.
        return false;
    }

    private static boolean isCandidateRecordType(SAMRecord record)
    {
        // Only care about primaries because only they have the read sequence.
        // And REDUX duplicate fragments don't provide additional support.
        return !record.isSecondaryOrSupplementary() && !record.getDuplicateReadFlag();
    }

    private boolean hasCandidateSoftClip(SAMRecord record)
    {
        boolean hasSignificantClip = leftSoftClipLength(record) >= mMinSoftClipBases || rightSoftClipLength(record) >= mMinSoftClipBases;
        return hasSignificantClip && clippedBasesAreCandidate(record);
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

    private static boolean isReduxUnmapped(SAMRecord record)
    {
        return record.hasAttribute(UNMAP_ATTRIBUTE);
    }

    @Nullable
    private static String reduxUnmappedOriginalContig(SAMRecord record)
    {
        String unmapCoords = record.getStringAttribute(UNMAP_ATTRIBUTE);
        if(unmapCoords == null)
        {
            return null;
        }
        int delimiterIndex = unmapCoords.lastIndexOf(UNMAPP_COORDS_DELIM);
        if(delimiterIndex <= 0)
        {
            throw new UserInputError("Malformed REDUX " + UNMAP_ATTRIBUTE + " attribute: " + unmapCoords);
        }
        return unmapCoords.substring(0, delimiterIndex);
    }

    private boolean isReduxUnmappedFromVirusDecoy(SAMRecord record)
    {
        if(!isReduxUnmapped(record))
        {
            return false;
        }
        String originalContig = reduxUnmappedOriginalContig(record);
        return originalContig != null && isVirusDecoyContig(originalContig);
    }

    private boolean isVirusDecoyContig(String contig)
    {
        return mRefVirusContigs.contains(contig);
    }

    private boolean hasViralKmer(SAMRecord record)
    {
        return mKmerIndex == null || mKmerIndex.hasViralKmer(record.getReadBases());
    }
}
