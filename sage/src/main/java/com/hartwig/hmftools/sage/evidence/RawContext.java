package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.SageConstants.MAX_SOFT_CLIP_LOW_QUAL_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SOFT_CLIP_HIGH_QUAL_PERC;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SOFT_CLIP_MIN_BASE_QUAL;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.ALIGNED;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.DELETED;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.LOW_QUAL;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.SKIPPED;
import static com.hartwig.hmftools.sage.evidence.VariantReadPositionType.SOFT_CLIP;

import com.hartwig.hmftools.sage.common.SimpleVariant;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class RawContext
{
    public final int ReadVariantIndex;
    public final VariantReadPositionType PositionType;

    protected static final RawContext INVALID_CONTEXT = new RawContext(-1, VariantReadPositionType.NONE);

    public RawContext(final int readVariantIndex, final VariantReadPositionType positionType)
    {
        ReadVariantIndex = readVariantIndex;
        PositionType = positionType;
    }

    public String toString()
    {
        return format("readVarIndex(%d) posType(%s)", ReadVariantIndex, PositionType);
    }

    public static RawContext createFromRead(final SimpleVariant variant, final SAMRecord record)
    {
        final Cigar cigar = record.getCigar();

        int readIndex = 0;
        int refBase = record.getAlignmentStart();
        RawContext rawContext = null;

        for(int i = 0; i < cigar.numCigarElements(); i++)
        {
            final CigarElement element = cigar.getCigarElement(i);

            switch(element.getOperator())
            {
                case H: // ignore hard clips - no need to skip either bases or positions
                case P: // ignore pads
                    break;

                case S:
                    if(i == 0)
                    {
                        rawContext = handleLeftSoftClip(variant, record, element);
                    }
                    else if(i == cigar.numCigarElements() - 1)
                    {
                        rawContext = handleRightSoftClip(variant, record, element, refBase);
                    }
                    readIndex += element.getLength();
                    break;

                case N:
                    rawContext = handleSkippedReference(variant, record, element, readIndex - 1, refBase - 1);
                    refBase += element.getLength();
                    break;

                case D:
                    rawContext = handleDelete(variant, element, readIndex - 1, refBase - 1);
                    refBase += element.getLength();
                    break;

                case I:

                    if(readIndex == 0) // unusual and ref base is still assumed to be the base prior
                        rawContext = handleInsert(variant, 0, refBase - 1);
                    else
                        rawContext = handleInsert(variant, readIndex - 1, refBase - 1);

                    readIndex += element.getLength();
                    break;

                case M:
                case EQ:
                case X:
                    rawContext = handleAlignment(variant, element, readIndex, refBase);
                    readIndex += element.getLength();
                    refBase += element.getLength();
                    break;

                default:
                    throw new IllegalStateException("Case statement didn't deal with op: " + element.getOperator() + "in CIGAR: " + cigar);
            }

            if(rawContext != null)
                break;
        }

        return rawContext != null ? rawContext : INVALID_CONTEXT;
    }

    private static RawContext handleLeftSoftClip(final SimpleVariant variant, final SAMRecord record, final CigarElement element)
    {
        if(variant.Position >= record.getAlignmentStart())
            return null;

        // set read index assuming a REF match, not a variant match
        // eg if soft-clip length = 10 (index 0-9), variant position is 96 vs start of 100, pos diff = 4, then read index = 10 - 4 = 6
        int variantPosDiff = record.getAlignmentStart() - variant.Position;
        int readIndex = element.getLength() - variantPosDiff;

        if(exceedsSoftClipLowBaseQual(record.getBaseQualities(), 0, element.getLength()))
            return new RawContext(readIndex, LOW_QUAL);

        return new RawContext(readIndex, SOFT_CLIP);
    }

    private static RawContext handleRightSoftClip(
            final SimpleVariant variant, final SAMRecord record, final CigarElement element, int refPosition)
    {
        int unclippedEnd = refPosition + element.getLength() - 1;

        if(variant.Position > unclippedEnd)
            return null;

        int scStartIndex = record.getBaseQualities().length - element.getLength();

        int variantPosDiff = variant.Position - record.getAlignmentEnd();
        int variantReadIndex = scStartIndex + variantPosDiff - 1;

        if(exceedsSoftClipLowBaseQual(record.getBaseQualities(), scStartIndex, element.getLength()))
            return new RawContext(variantReadIndex, LOW_QUAL);

        return new RawContext(variantReadIndex, SOFT_CLIP);
    }

    private static RawContext handleAlignment(final SimpleVariant variant, final CigarElement element, int readIndex, int refPosition)
    {
        // when the variant is an indel, its position will be at the last aligned base of the previous M element
        // in this case these need to be ignored for the purpose of establishing raw support for the variant
        int refPositionEnd = refPosition + element.getLength() - 1;

        if(variant.isIndel())
            --refPositionEnd;

        if(refPosition <= variant.Position && variant.Position <= refPositionEnd)
        {
            int readIndexOffset = variant.Position - refPosition;
            int variantReadIndex = readIndex + readIndexOffset;

            return new RawContext(variantReadIndex, ALIGNED);
        }

        return null;
    }

    private static RawContext handleInsert(final SimpleVariant variant, int readIndex, int refPosition)
    {
        return refPosition == variant.Position ? new RawContext(readIndex, ALIGNED) : null;
    }

    private static RawContext handleDelete(final SimpleVariant variant, final CigarElement element, int readIndex, int refPosition)
    {
        int refPositionEnd = refPosition + element.getLength();
        if(refPosition == variant.Position)
        {
            return new RawContext(readIndex, ALIGNED);
        }
        else if(positionWithin(variant.Position, refPosition, refPositionEnd))
        {
            return new RawContext(readIndex, DELETED);
        }

        return null;
    }

    private static RawContext handleSkippedReference(
            final SimpleVariant variant, final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
    {
        int refPositionEnd = refPosition + element.getLength();
        if(refPositionEnd >= variant.Position)
        {
            return new RawContext(readIndex, SKIPPED);
        }

        return null;
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
