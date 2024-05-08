package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.SageConstants.MAX_SOFT_CLIP_LOW_QUAL_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SOFT_CLIP_HIGH_QUAL_PERC;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SOFT_CLIP_MIN_BASE_QUAL;
import static com.hartwig.hmftools.sage.candidate.RefContextConsumer.ignoreSoftClipAdapter;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.ALIGNED;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.DELETED;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.SKIPPED;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.SOFT_CLIP;

import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.SplitReadUtils;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class RawContextCigarHandler implements CigarHandler
{
    private final SimpleVariant mVariant;
    private final boolean mIsInsert;
    private final boolean mIsDelete;
    private final boolean mIsSNV;

    private RawContext mResult;

    public RawContextCigarHandler(final SimpleVariant variant)
    {
        mVariant = variant;
        mIsInsert = variant.isInsert();
        mIsDelete = variant.isDelete();
        mIsSNV = variant.ref().length() == variant.alt().length();
    }

    public RawContext result()
    {
        return mResult;
    }

    @Override
    public void handleLeftSoftClip(final SAMRecord record, final CigarElement element)
    {
        if(mVariant.position() >= record.getAlignmentStart())
            return;

        if(ignoreSoftClipAdapter(record))
            return;

        if(exceedsSoftClipLowBaseQual(record.getBaseQualities(), 0, element.getLength()))
            return;

        // set read index assuming a REF match, not a variant match
        // eg if soft-clip length = 10 (index 0-9), variant position is 96 vs start of 100, pos diff = 4, then read index = 10 - 4 = 6
        int variantPosDiff = record.getAlignmentStart() - mVariant.Position;
        int readIndex = element.getLength() - variantPosDiff;

        mResult = new RawContext(readIndex, SOFT_CLIP);
    }

    @Override
    public void handleRightSoftClip(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
    {
        if(mResult != null && !mIsInsert)
            return;

        int unclippedEnd = refPosition + element.getLength() - 1;

        if(mVariant.position() > unclippedEnd)
            return;

         if(ignoreSoftClipAdapter(record))
             return;

        int scStartIndex = record.getBaseQualities().length - element.getLength();

        if(exceedsSoftClipLowBaseQual(record.getBaseQualities(), scStartIndex, element.getLength()))
            return;

        int variantPosDiff = mVariant.Position - record.getAlignmentEnd();
        int variantReadIndex = scStartIndex + variantPosDiff - 1;
        mResult = new RawContext(variantReadIndex, SOFT_CLIP);

        /* CLEANUP: cannot explain any of the logic below

        if(mIsInsert)
        {
            // for inserts from the soft-clipped bases, the variant's position will be the last ref position prior to the SC start
            if(mVariant.position() != record.getAlignmentEnd())
                return;

            int readVariantStartPos = readIndex - 1;

            boolean altSupport = element.getLength() >= mVariant.alt().length() && matchesString(record, readVariantStartPos, mVariant.alt());

            if(!altSupport)
                return;

            mResult = new RawContext(readVariantStartPos, SOFT_CLIP);
        }
        else
        {
            if(mVariant.position() >= refPosition && mVariant.position() <= refPositionEnd)
            {
                int alignmentEnd = record.getAlignmentEnd();
                int actualIndex = record.getReadPositionAtReferencePosition(alignmentEnd) - 1 - alignmentEnd + mVariant.position();

                mResult = new RawContext(actualIndex, SOFT_CLIP);
            }
        }
        */
    }

    @Override
    public void handleAlignment(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
    {
        if(mResult != null)
            return;

        // when the variant is an indel, its position will be at the last aligned base of the previous M element
        // in this case these need to be ignored for the purpose of establishing raw support for the variant
        int refPositionEnd = refPosition + element.getLength() - 1;

        if(mVariant.isIndel())
            --refPositionEnd;

        if(refPosition <= mVariant.position() && mVariant.position() <= refPositionEnd)
        {
            int readIndexOffset = mVariant.position() - refPosition;
            int variantReadIndex = readIndex + readIndexOffset;

            mResult = new RawContext(variantReadIndex, ALIGNED);
        }
    }

    @Override
    public void handleInsert(final SAMRecord record, final CigarElement e, int readIndex, int refPosition)
    {
        if(mResult != null)
            return;

        if(refPosition == mVariant.position())
        {
            mResult = new RawContext(readIndex, ALIGNED);
        }
    }

    @Override
    public void handleDelete(final SAMRecord record, final CigarElement e, int readIndex, int refPosition)
    {
        if(mResult != null)
            return;

        int refPositionEnd = refPosition + e.getLength();
        if(refPosition == mVariant.position())
        {
            mResult = new RawContext(readIndex, ALIGNED);
        }
        else if(positionWithin(mVariant.position(), refPosition, refPositionEnd))
        {
            mResult = new RawContext(readIndex, DELETED);
        }
    }

    @Override
    public void handleSkippedReference(final SAMRecord record, final CigarElement e, int readIndex, int refPosition)
    {
        if(mResult != null)
            return;

        if(e.getLength() > SplitReadUtils.MAX_SKIPPED_REFERENCE_REGIONS)
        {
            int refPositionEnd = refPosition + e.getLength();
            if(refPositionEnd >= mVariant.position())
            {
                mResult = new RawContext(readIndex, SKIPPED);
            }
        }

        handleDelete(record, e, readIndex, refPosition);
    }

    private static boolean matchesString(final SAMRecord record, int index, final String expected)
    {
        if(index < 0)
            return false;

        int readLength = record.getReadBases().length;

        for(int i = 0; i < expected.length(); i++)
        {
            if(index + i >= readLength)
                return false;

            if((byte) expected.charAt(i) != record.getReadBases()[index + i])
                return false;
        }

        return true;
    }

    public static boolean exceedsSoftClipLowBaseQual(final byte[] baseQualities, int startIndex, int scLength)
    {
        if(scLength < MAX_SOFT_CLIP_LOW_QUAL_COUNT)
            return false;

        double requiredHighQual = MIN_SOFT_CLIP_HIGH_QUAL_PERC * scLength;
        int aboveQual = 0;

        for(int i = startIndex; i < startIndex + scLength; ++i)
        {
            if(baseQualities[i] >= MIN_SOFT_CLIP_MIN_BASE_QUAL)
            {
                ++aboveQual;

                if(aboveQual >= requiredHighQual)
                    return false;
            }
        }

        int lowQualCount = scLength - aboveQual;
        return lowQualCount >= MAX_SOFT_CLIP_LOW_QUAL_COUNT;
    }
}
