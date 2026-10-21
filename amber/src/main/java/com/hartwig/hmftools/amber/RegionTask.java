package com.hartwig.hmftools.amber;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;

import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class RegionTask
{
    private final PositionEvidenceChecker mEvidenceChecker;
    public final ChrBaseRegion Region;

    private final List<PositionEvidence> mPositions;
    private int mCurrentIndex;
    private boolean mComplete;

    public RegionTask(final PositionEvidenceChecker evidenceChecker, final String chromosome, final PositionEvidence baseDepth)
    {
        mEvidenceChecker = evidenceChecker;
        Region = new ChrBaseRegion(chromosome, baseDepth.Position, baseDepth.Position);
        mPositions = Lists.newArrayList(baseDepth);
        mCurrentIndex = 0;
        mComplete = false;
    }

    public void addPosition(final PositionEvidence posEvidence)
    {
        mPositions.add(posEvidence);
        Region.setEnd(max(Region.end(), posEvidence.Position));
    }

    public void processRecord(final SAMRecord record)
    {
        int[] fragmentBoundaries = checkFragmentBoundaries(record);

        int alignmentStart, alignmentEnd;
        int readAlignmentStart = record.getAlignmentStart();

        if(fragmentBoundaries != null)
        {
            alignmentStart = fragmentBoundaries[0];
            alignmentEnd = fragmentBoundaries[1];

            if(alignmentStart > alignmentEnd)
                return;
        }
        else
        {
            alignmentStart = readAlignmentStart;
            alignmentEnd = record.getAlignmentEnd();
        }

        int index = mCurrentIndex;

        ReadCigarState cigarState = null;

        for(; index < mPositions.size(); ++index)
        {
            PositionEvidence posEvidence = mPositions.get(index);

            if(readAlignmentStart > posEvidence.Position) // permanently move past these positions based on BAM alignment
            {
                ++mCurrentIndex;
                continue;
            }

            if(alignmentStart > posEvidence.Position) // this fragment's upper read may not process this position, but other reads could
                continue;

            if(alignmentEnd < posEvidence.Position)
                break;

            if(cigarState == null)
            {
                // use cigar to move through read bases efficiently
                cigarState = new ReadCigarState(record);
            }

            cigarState.moveToPosition(posEvidence.Position);

            if(!cigarState.isValid())
                break;

            mEvidenceChecker.addEvidence(posEvidence, record, cigarState.ReadIndex, cigarState.isIndel());
        }

        if(mCurrentIndex >= mPositions.size())
            mComplete = true;
    }

    public boolean isComplete() { return mComplete; }

    public int positionCount() { return mPositions.size(); }

    public String toString()
    {
        return format("region(%s) positions(%d) index(%d)", Region, mPositions.size(), mCurrentIndex);
    }

    @VisibleForTesting
    public static int[] checkFragmentBoundaries(final SAMRecord record)
    {
        // applies to concordant overlapping reads only
        if(!record.getReadPairedFlag() || record.getMateUnmappedFlag())
            return null;

        if(record.getReferenceIndex() != record.getMateReferenceIndex())
            return null;

        if(record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
            return null;

        if(record.getAlignmentEnd() < record.getMateAlignmentStart())
            return null;

        String mateCigar = record.getStringAttribute(MATE_CIGAR_ATTRIBUTE);

        if(mateCigar == null)
            return null;

        int readStart = record.getAlignmentStart();
        int readEnd = record.getAlignmentEnd();

        int mateStart = record.getMateAlignmentStart();

        List<CigarElement> mateCigarElements = cigarFromStr(mateCigar).getCigarElements();

        int mateEnd = mateStart + mateCigarElements.stream()
                .filter(x -> x.getOperator().consumesReferenceBases()).mapToInt(x -> x.getLength()).sum() - 1;

        if(mateEnd < readStart)
            return null;

        boolean isFirstRead;

        if(readStart < mateStart)
        {
            isFirstRead = true;
        }
        else if(readStart == mateStart)
        {
            if(readEnd == mateEnd)
                isFirstRead = !record.getReadNegativeStrandFlag(); // as per BAM coordinate sorting
            else
                isFirstRead = readEnd < mateEnd; // so as to process the shorter read first
        }
        else
        {
            isFirstRead = false;
        }

        // avoid processing the same bases from the earlier read twice
        if(!isFirstRead)
        {
            readStart = max(readStart, mateEnd + 1);
        }

        // limit the 3' end of reads to the unclipped part of the mate's 5' end to ignore adapter bases
        if(record.getReadNegativeStrandFlag())
        {
            int unclippedMateStart = mateStart;

            if(mateCigarElements.get(0).getOperator() == S)
                unclippedMateStart -= mateCigarElements.get(0).getLength();

            readStart = max(unclippedMateStart, readStart);
        }
        else
        {
            int unclippedMateEnd = mateEnd;

            int lastIndex = mateCigarElements.size() - 1;
            if(mateCigarElements.get(lastIndex).getOperator() == S)
                unclippedMateEnd += mateCigarElements.get(lastIndex).getLength();

            readEnd = min(unclippedMateEnd, readEnd);
        }

        return new int[] { readStart, readEnd };
    }
}
