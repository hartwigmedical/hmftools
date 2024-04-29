package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarStringFromElements;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.sage.evidence.ArtefactContext;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class VariantReadContextBuilder
{
    private final int mFlankSize;

    public VariantReadContextBuilder(int flankSize)
    {
        mFlankSize = flankSize;
    }

    public VariantReadContext createContext(
            final SimpleVariant variant, final SAMRecord read, int varReadIndex, final RefSequence refSequence)
    {
        try
        {
            VariantReadContext readContext = buildContext(variant, read, varReadIndex, refSequence);

            if(readContext == null)
                return null;

            // enforce full flanks
            if(readContext.leftFlankLength() < mFlankSize || readContext.rightFlankLength() < mFlankSize)
                return null;

            readContext.setArtefactContext(ArtefactContext.buildContext(readContext));

            return readContext;
        }
        catch(Exception e)
        {
            SG_LOGGER.error("var({}) error building readContext, varReadIndex({}) read({}) error: {}",
                    variant, varReadIndex, readToString(read), e.toString());
            e.printStackTrace();

            return null;
        }
    }

    private VariantReadContext buildContext(
            final SimpleVariant variant, final SAMRecord read, int varIndexInRead, final RefSequence refSequence)
    {
        /* Routine:
            - start with a basic core around the variant
            - test for homology and right-align if required
            - continue in each direction until repeat ends
            - add core bases then flank
            - extra ref bases to cover the core length
        */

        int readCoreStart, readCoreEnd;
        Microhomology homology = null;

        if(variant.isIndel())
        {
            readCoreStart = varIndexInRead - MIN_CORE_DISTANCE + 1;
            readCoreEnd = varIndexInRead + (variant.isInsert() ? variant.indelLength() + 1 : 1) + MIN_CORE_DISTANCE - 1;

            homology = findHomology(variant, read, varIndexInRead);

            if(homology != null)
                readCoreEnd += homology.Length;
        }
        else
        {
            readCoreStart = varIndexInRead - MIN_CORE_DISTANCE;
            readCoreEnd = varIndexInRead + variant.Alt.length() - 1 + MIN_CORE_DISTANCE;
        }

        if(readCoreStart < 0 || readCoreEnd >= read.getReadBases().length)
            return null;

        final byte[] readBases = read.getReadBases();

        RepeatBoundaries repeatBoundaries = findRepeatBoundaries(readCoreStart, readCoreEnd, readBases);

        RepeatInfo maxRepeat = null;

        if(repeatBoundaries != null)
        {
            readCoreStart = min(readCoreStart, repeatBoundaries.LowerIndex);
            readCoreEnd = max(readCoreEnd, repeatBoundaries.UpperIndex);
            maxRepeat = repeatBoundaries.MaxRepeat;
        }

        int readFlankStart = max(readCoreStart - mFlankSize, 0);
        int readFlankEnd = min(readCoreEnd + mFlankSize, readBases.length - 1);

        if(readFlankStart < 0 || readFlankEnd >= read.getReadBases().length)
            return null;

        // build a CIGAR from the read to cover this range
        ReadCigarInfo readCigarInfo = buildReadCigar(read, readFlankStart, readFlankEnd);

        readFlankStart = readCigarInfo.FlankIndexStart;
        readFlankEnd = readCigarInfo.FlankIndexEnd;

        byte[] contextReadBases = Arrays.subsetArray(readBases, readFlankStart, readFlankEnd);

        int readContextOffset = readFlankStart;
        int coreIndexStart = readCoreStart - readContextOffset;
        int readVarIndex = varIndexInRead - readContextOffset;
        int coreIndexEnd = readCoreEnd - readContextOffset;

        int alignmentStart = max(read.getAlignmentStart(), readCigarInfo.UnclippedStart);
        int alignmentEnd = min(read.getAlignmentEnd(), readCigarInfo.UnclippedEnd);

        // ref bases are the core width around the variant's position
        int leftCoreLength = readVarIndex - coreIndexStart;
        int coreLength = coreIndexEnd - coreIndexStart + 1;
        int refPosStart = variant.Position - leftCoreLength;
        int refPosEnd;

        if(variant.isDelete())
        {
            // ensure is long enough to cover the ref bases prior to deletion, so find the core position end in the ref
            int rightCoreLength = coreIndexEnd - readVarIndex;
            int rightCoreLengthExtension = max(rightCoreLength - MIN_CORE_DISTANCE, 0);
            refPosEnd = variant.positionEnd() + 1 + rightCoreLengthExtension;
        }
        else
        {
            refPosEnd = refPosStart + coreLength - 1;
        }

        byte[] refBases = refSequence.baseRange(refPosStart, refPosEnd);

        return new VariantReadContext(
                variant, alignmentStart, alignmentEnd, refBases, contextReadBases, readCigarInfo.Cigar,
                coreIndexStart, readVarIndex, coreIndexEnd, homology, maxRepeat);
    }

    private RepeatBoundaries findRepeatBoundaries(int readCoreStart, int readCoreEnd, final byte[] readBases)
    {
        return RepeatBoundaries.findRepeatBoundaries(readBases, readCoreStart, readCoreEnd, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);
    }

    private class ReadCigarInfo
    {
        public List<CigarElement> Cigar;
        public final int UnclippedStart;  // note: these are unclipped values
        public final int UnclippedEnd;

        // stored for the scenario where an indel in the flanks pushes out the alignment beyond the standard flank length
        public final int FlankIndexStart;
        public final int FlankIndexEnd;

        public ReadCigarInfo(
                final List<CigarElement> cigar, final int unclippedStart, final int unclippedEnd,
                final int flankIndexStart, final int flankIndexEnd)
        {
            Cigar = cigar;
            UnclippedStart = unclippedStart;
            UnclippedEnd = unclippedEnd;
            FlankIndexStart = flankIndexStart;
            FlankIndexEnd = flankIndexEnd;
        }

        public String toString()
        {
            return format("%s align(%d-%d) flankIndex(%d-%d)",
                    cigarStringFromElements(Cigar), UnclippedStart, UnclippedEnd, FlankIndexStart, FlankIndexEnd);
        }
    }

    private ReadCigarInfo buildReadCigar(final SAMRecord read, int indexStart, int indexEnd)
    {
        // builds a cigar around the specified index boundaries and calculates the corresponding alignment positions
        List<CigarElement> cigar = Lists.newArrayList();

        int readIndex = 0;
        int refPosition = read.getAlignmentStart();
        int unclippedPosStart = 0;
        int unclippedPosEnd = 0;

        int finalIndexStart = indexStart;
        int finalIndexEnd = indexEnd;

        for(CigarElement element : read.getCigar().getCigarElements())
        {
            if(readIndex == 0 && element.getOperator() == S)
            {
                // set to unclipped ref position so the alignment can capture corresponding ref bases
                refPosition -= element.getLength();
            }

            int elementEndIndex = readIndex + element.getLength() - 1;

            if(elementEndIndex >= indexStart)
            {
                int elementStart = max(readIndex, indexStart);
                int elementEnd = min(elementEndIndex, indexEnd);

                cigar.add(new CigarElement(elementEnd - elementStart + 1, element.getOperator()));

                if(unclippedPosStart == 0)
                {
                    if(element.getOperator() == I)
                    {
                        // handles an insert that pushes the alignment out - always take the prior alignment base and reduce index start
                        // eg looking to find the alignment boundary for index start 12 for 10M5I1350M, alignment start == 100
                        // so at the insert element, read index = 10, ref pos = 110 (pointing at next ref base)
                        int extraIndexStart = indexStart - readIndex + 1;
                        cigar.add(0, new CigarElement(1, M));

                        finalIndexStart -= extraIndexStart;
                        unclippedPosStart = max(refPosition - 1, read.getAlignmentStart());
                    }
                    else if(element.getOperator() == D)
                    {
                        unclippedPosStart = max(refPosition - 1, read.getAlignmentStart());
                    }
                    else
                    {
                        unclippedPosStart = refPosition + (indexStart - readIndex);
                    }
                }

                if(unclippedPosEnd == 0 && elementEndIndex >= indexEnd)
                {
                    if(element.getOperator() == I)
                    {
                        // similar extension to the above
                        int extraIndexEnd = elementEndIndex + 1 - indexEnd;
                        cigar.add(new CigarElement(1, M));

                        finalIndexEnd += extraIndexEnd;

                        unclippedPosEnd = refPosition; // already pointing at the next M (aligned) base
                    }
                    else if(element.getOperator() == D)
                    {
                        unclippedPosEnd = refPosition + element.getLength();
                    }
                    else
                    {
                        unclippedPosEnd = refPosition + (indexEnd - readIndex); // unclipped so not bound by the read's alignment
                    }
                }
            }

            if(element.getOperator().consumesReadBases())
                readIndex += element.getLength();

            if(element.getOperator().consumesReferenceBases() || element.getOperator() == S)
                refPosition += element.getLength();

            if(readIndex > indexEnd)
                break;
        }

        return new ReadCigarInfo(cigar, unclippedPosStart, unclippedPosEnd, finalIndexStart, finalIndexEnd);
    }

    public static Microhomology findHomology(final SimpleVariant variant, final SAMRecord read, int varReadIndex)
    {
        if(!variant.isIndel())
            return null;

        int indelAltLength = abs(variant.indelLength());

        StringBuilder homology = null;
        String indelBases = variant.isInsert() ? variant.alt().substring(1) : variant.ref().substring(1);

        final byte[] readBases = read.getReadBases();

        // start looking in the read in the first base after the variant

        int homReadIndexStart = variant.isInsert() ? varReadIndex + indelAltLength + 1 : varReadIndex + 1;
        int homReadIndex = homReadIndexStart;

        for(int i = 0; i < indelBases.length() && homReadIndex < readBases.length; ++i, ++homReadIndex)
        {
            byte indelBase = (byte)indelBases.charAt(i);
            byte postIndelReadBase = readBases[homReadIndex];

            if(indelBase != postIndelReadBase)
                break;

            if(homology == null)
                homology = new StringBuilder();

            homology.append((char)indelBase);
        }

        if(homology == null)
            return null;

        String homologyBases = homology.toString();
        int homologyLength = homologyBases.length();

        if(homologyLength == indelAltLength && homReadIndex < readBases.length)
        {
            // continue searching for repeats of the homology to find the point of right-alignment, allow partials at the end
            boolean matched = true;

            while(matched)
            {
                int matchCount = 0;

                for(int i = 0; i < indelAltLength && homReadIndex < readBases.length; ++i, ++homReadIndex)
                {
                    byte homologyBase = (byte)homologyBases.charAt(i);
                    byte postIndelReadBase = readBases[homReadIndex];

                    if(homologyBase != postIndelReadBase)
                    {
                        matched = false;
                        break;
                    }

                    ++matchCount;
                }

                if(matchCount > 0)
                    homologyLength += matchCount;
                else
                    break;
            }
        }

        return new Microhomology(homologyBases, homologyLength);
    }

    public static final int INVALID_INDEX_POS = -1;

    public static int findPositionStart(
            int variantPosition, int leftCoreLength, int alignmentStart, final List<CigarElement> readCigar, int readIndex)
    {
        if(readCigar.size() == 1)
            return variantPosition - leftCoreLength;

        int position = findReadPositionFromIndex(alignmentStart, readCigar, readIndex);

        return position > 0 ? position : variantPosition - leftCoreLength;
    }

    public static int findPositionEnd(
            int variantPosition, int rightCoreLength, int alignmentStart, final List<CigarElement> readCigar, int readIndex)
    {
        if(readCigar.size() == 1)
            return variantPosition + rightCoreLength;

        int position = findReadPositionFromIndex(alignmentStart, readCigar, readIndex);

        return position > 0 ? position : variantPosition + rightCoreLength;
    }

    public static int findReadPositionFromIndex(int alignmentStart, final List<CigarElement> readCigar, int readIndex)
    {
        int refPosition = alignmentStart;
        int index = 0;

        for(CigarElement element : readCigar)
        {
            if(index + element.getLength() >= readIndex && element.getOperator().consumesReadBases())
            {
                if(element.getOperator().consumesReferenceBases())
                    refPosition += readIndex - index;

                return refPosition;
            }

            if(element.getOperator().consumesReferenceBases())
                refPosition += element.getLength();

            if(element.getOperator().consumesReadBases())
                index += element.getLength();
        }

        // shouldn't occur
        return INVALID_INDEX_POS;
    }

    public static int determineUpperAltIndex(
            final SimpleVariant variant, final byte[] readBases, final byte[] refBases,
            final int varReadIndex, final int varRefIndex, final int coreIndexEnd)
    {
        // find the first base of difference (ref vs alt) up and down from the variant's position, and cap at the core indices
        if(!variant.isIndel())
        {
            return varReadIndex + variant.Alt.length() - 1;
        }

        int refIndex = varRefIndex;
        int readIndex = varReadIndex;

        for(; readIndex <= coreIndexEnd & refIndex < refBases.length; ++readIndex, ++refIndex)
        {
            if(refBases[refIndex] != readBases[readIndex])
                break;
        }

        int upperRefIndex = variant.isInsert() ? varReadIndex + variant.Alt.length() : varReadIndex + 1;

        return max(min(readIndex, coreIndexEnd), upperRefIndex);
    }

}
