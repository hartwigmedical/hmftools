package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.MAX_SOFT_CLIP_LOW_QUAL_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SOFT_CLIP_HIGH_QUAL_PERC;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SOFT_CLIP_MIN_BASE_QUAL;
import static com.hartwig.hmftools.sage.candidate.RefContextConsumer.ignoreSoftClipAdapter;

import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.read.SplitReadUtils;

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
        if(mVariant.position() > record.getAlignmentStart())
            return;

        if(ignoreSoftClipAdapter(record))
            return;

        if(exceedsSoftClipLowBaseQual(record.getBaseQualities(), 0, element.getLength()))
            return;

        int readStartPos = record.getReadPositionAtReferencePosition(record.getAlignmentStart());
        int readIndex = readStartPos - 1 - record.getAlignmentStart()
                + mVariant.position() - mVariant.alt().length() + mVariant.ref().length();

        boolean altSupport = mIsInsert && element.getLength() >= mVariant.alt().length() && matchesString(record, readIndex, mVariant.alt());
        int baseQuality = altSupport ? avgBaseQuality(readIndex, record, mVariant.alt().length()) : 0;

        mResult = RawContext.inSoftClip(readIndex, altSupport, baseQuality);
    }

    @Override
    public void handleRightSoftClip(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
    {
        if(mResult != null && !mIsInsert)
            return;

        int refPositionEnd = refPosition + element.getLength() - 1;
        if(refPositionEnd < mVariant.position())
            return;

         if(ignoreSoftClipAdapter(record))
             return;

        int scStartIndex = record.getBaseQualities().length - element.getLength();

        if(exceedsSoftClipLowBaseQual(record.getBaseQualities(), scStartIndex, element.getLength()))
            return;

        if(mIsInsert)
        {
            // for inserts from the soft-clipped bases, the variant's position will be the last ref position prior to the SC start
            if(mVariant.position() != record.getAlignmentEnd())
                return;

            int readVariantStartPos = readIndex - 1;

            boolean altSupport = element.getLength() >= mVariant.alt().length() && matchesString(record, readVariantStartPos, mVariant.alt());
            int baseQuality = altSupport ? avgBaseQuality(readVariantStartPos, record, mVariant.alt().length()) : 0;

            if(!altSupport)
                return;

            mResult = RawContext.inSoftClip(readVariantStartPos, altSupport, baseQuality);
        }
        else
        {
            if(mVariant.position() >= refPosition && mVariant.position() <= refPositionEnd)
            {
                int alignmentEnd = record.getAlignmentEnd();
                int actualIndex = record.getReadPositionAtReferencePosition(alignmentEnd) - 1 - alignmentEnd + (int) mVariant.position();

                mResult = RawContext.inSoftClip(actualIndex, false, 0);
            }
        }
    }

    @Override
    public void handleAlignment(final SAMRecord record, final CigarElement element, boolean beforeIndel, int readIndex, int refPosition)
    {
        if(mResult != null)
            return;

        // check carefully when element is followed by an INDEL
        int refPositionEnd = refPosition + element.getLength() - 1;

        if(refPosition <= mVariant.position() && mVariant.position() <= refPositionEnd)
        {
            int readIndexOffset = mVariant.position() - refPosition;
            int variantReadIndex = readIndex + readIndexOffset;

            int baseQuality = record.getBaseQualities()[variantReadIndex];

            // this logic is inadequate for INDELs but is handled in the index matching routine instead - ie to reverse any incorrect
            // ref support here
            boolean altSupport = mIsSNV && refPositionEnd >= mVariant.end() && matchesString(record, variantReadIndex, mVariant.alt());
            boolean refSupport = !altSupport && matchesFirstBase(record, variantReadIndex, mVariant.ref(), true);
            mResult = RawContext.alignment(variantReadIndex, altSupport, refSupport, baseQuality);
        }
    }

    @Override
    public void handleInsert(final SAMRecord record, final CigarElement e, int readIndex, int refPosition)
    {
        if(mResult != null)
            return;

        if(refPosition == mVariant.position())
        {
            boolean altSupport = mIsInsert && e.getLength() == mVariant.alt().length() - 1 && matchesString(record, readIndex, mVariant.alt());
            int baseQuality = altSupport ? avgBaseQuality(readIndex, record, mVariant.alt().length()) : 0;
            mResult = RawContext.indel(readIndex, altSupport, baseQuality);
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
            boolean altSupport = mIsDelete && e.getLength() == mVariant.ref().length() - 1
                    && matchesFirstBase(record, readIndex, mVariant.ref(), false);

            int baseQuality = altSupport ? avgBaseQuality(readIndex, record, 2) : 0;
            mResult = RawContext.indel(readIndex, altSupport, baseQuality);
        }
        else if(positionWithin(mVariant.position(), refPosition, refPositionEnd))
        {
            mResult = RawContext.inDelete(readIndex);
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
                mResult = RawContext.inSkipped(readIndex);
            }
        }

        handleDelete(record, e, readIndex, refPosition);
    }

    private static boolean matchesFirstBase(final SAMRecord record, int index, final String expected, boolean allowLowQual)
    {
        if(expected.charAt(0) == record.getReadBases()[index])
            return true;

        return allowLowQual && record.getBaseQualities()[index] < MATCHING_BASE_QUALITY;
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

    private int avgBaseQuality(int readIndex, SAMRecord record, int length)
    {
        int maxIndex = Math.min(readIndex + length, record.getBaseQualities().length) - 1;

        double qualityTotal = 0;

        for(int i = readIndex; i <= maxIndex; i++)
        {
            qualityTotal += record.getBaseQualities()[i];
        }

        int baseLength = maxIndex - readIndex + 1;
        return qualityTotal > 0 ? (int)(qualityTotal / baseLength) : 0;
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
