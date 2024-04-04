package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.sage.old.ReadContext;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class VariantReadContextBuilder
{
    private final int mFlankSize;

    public VariantReadContextBuilder(int flankSize)
    {
        mFlankSize = flankSize;
    }

    private static final int REPEAT_SEARCH_LENGTH = MAX_REPEAT_LENGTH * 2;

    private class RepeatBoundaryInfo
    {
        public final int CoreStart;
        public final int CoreEnd;
        public final RepeatInfo MaxRepeat;

        public RepeatBoundaryInfo(final int coreStart, final int coreEnd, final RepeatInfo maxRepeat)
        {
            CoreStart = coreStart;
            CoreEnd = coreEnd;
            MaxRepeat = maxRepeat;
        }
    }

    private RepeatBoundaryInfo findRepeatBoundaries(int lowerRefIndex, int upperRefIndex, final byte[] readBases)
    {
        int repeatSearchStart = lowerRefIndex - REPEAT_SEARCH_LENGTH;

        RepeatInfo lowerRepeat = RepeatInfo.findMaxRepeat(
                readBases, repeatSearchStart, readBases.length - 1,
                MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT, true, lowerRefIndex);

        int coreStart = lowerRefIndex;
        int coreEnd = upperRefIndex;

        RepeatInfo maxRepeat = lowerRepeat;

        if(lowerRepeat != null)
        {
            coreStart = lowerRepeat.Index - 1;
        }

        if(lowerRepeat == null || lowerRepeat.endIndex() <= upperRefIndex)
        {
            RepeatInfo upperRepeat = RepeatInfo.findMaxRepeat(
                    readBases, upperRefIndex, readBases.length - 1,
                    MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT, true, upperRefIndex);

            if(upperRepeat != null)
            {
                coreEnd = upperRepeat.endIndex() + 1;

                if(maxRepeat == null || upperRepeat.repeatLength() > maxRepeat.repeatLength())
                    maxRepeat = upperRepeat;
            }
        }

        return new RepeatBoundaryInfo(coreStart, coreEnd, maxRepeat);
    }

    private class ReadCigarInfo
    {
        public List<CigarElement> Cigar;
        public final int AlignmentStart;
        public final int AlignmentEnd;
        public final int FlankIndexStart;
        public final int FlankIndexEnd;

        public ReadCigarInfo(
                final List<CigarElement> cigar, final int alignmentStart, final int alignmentEnd,
                final int flankIndexStart, final int flankIndexEnd)
        {
            Cigar = cigar;
            AlignmentStart = alignmentStart;
            AlignmentEnd = alignmentEnd;
            FlankIndexStart = flankIndexStart;
            FlankIndexEnd = flankIndexEnd;
        }
    }

    private ReadCigarInfo buildReadCigar(final SAMRecord read, int indexStart, int indexEnd)
    {
        List<CigarElement> cigar = Lists.newArrayList();

        int readIndex = 0;
        int refPosition = read.getAlignmentStart();
        int alignmentStart = 0;
        int alignmentEnd = 0;

        for(CigarElement element : read.getCigar().getCigarElements())
        {
            int elementEndIndex = readIndex + element.getLength() - 1;

            if(elementEndIndex >= indexStart)
            {
                int elementStart = max(readIndex, indexStart);
                int elementEnd = min(elementEndIndex, indexEnd);

                cigar.add(new CigarElement(elementEnd - elementStart + 1, element.getOperator()));

                if(alignmentStart == 0)
                {
                    if(element.getOperator().isIndel())
                        alignmentStart = max(refPosition - 1, read.getAlignmentStart());
                    else
                        alignmentStart = refPosition + (indexStart - readIndex);
                }

                if(alignmentEnd == 0 && elementEndIndex >= indexEnd)
                {
                    if(element.getOperator() == I)
                        alignmentEnd = refPosition + 1;
                    if(element.getOperator() == D)
                        alignmentEnd = refPosition + element.getLength();
                    else
                        alignmentEnd = min(refPosition + (indexEnd - readIndex), read.getAlignmentEnd());
                }
            }

            if(element.getOperator().consumesReadBases())
                readIndex += element.getLength();

            if(element.getOperator().consumesReferenceBases())
                refPosition += element.getLength();

            if(readIndex > indexEnd)
                break;
        }

        return new ReadCigarInfo(cigar, alignmentStart, alignmentEnd, indexStart, indexEnd);
    }

    public VariantReadContext createMnvContext(
            final SimpleVariant variant, final SAMRecord read, int varReadIndex, final RefSequence refSequence)
    {
        int readCoreStart = varReadIndex - MIN_CORE_DISTANCE;
        int readCoreEnd = varReadIndex + variant.Alt.length() - 1 + MIN_CORE_DISTANCE;

        final byte[] readBases = read.getReadBases();

        RepeatBoundaryInfo repeatBoundaryInfo = findRepeatBoundaries(readCoreStart, readCoreEnd, readBases);

        readCoreStart = repeatBoundaryInfo.CoreStart;
        readCoreEnd = repeatBoundaryInfo.CoreEnd;

        int readFlankStart = max(readCoreStart - mFlankSize, 0);
        int readFlankEnd = min(readCoreEnd + mFlankSize, readBases.length - 1);

        // build a CIGAR from the read to cover this range
        ReadCigarInfo readCigarInfo = buildReadCigar(read, readFlankStart, readFlankEnd);

        readFlankStart = readCigarInfo.FlankIndexStart;
        readFlankEnd = readCigarInfo.FlankIndexEnd;

        byte[] contextReadBases = Arrays.subsetArray(readBases, readFlankStart, readFlankEnd);
        byte[] refBases = refSequence.baseRange(readCigarInfo.AlignmentStart, readCigarInfo.AlignmentEnd);

        int readContextOffset = readFlankStart;
        int coreIndexStart = readCoreStart - readContextOffset;
        int readVarIndex = varReadIndex - readContextOffset;
        int coreIndexEnd = readCoreEnd - readContextOffset;

        return new VariantReadContext(
                variant, readCigarInfo.AlignmentStart, readCigarInfo.AlignmentEnd, refBases,
                contextReadBases, readCigarInfo.Cigar, coreIndexStart, readVarIndex, coreIndexEnd, "", repeatBoundaryInfo.MaxRepeat);
    }

    public static Microhomology findHomology(final SimpleVariant variant, final SAMRecord read, int varReadIndex, final RefSequence refSequence)
    {
        if(!variant.isIndel())
            return null;

        int indelAltLength = abs(variant.indelLength());

        StringBuilder homology = null;
        final byte[] readBases = read.getReadBases();
        int refBaseStartPos = variant.isInsert() ? variant.Position + 1 : variant.Position + indelAltLength + 1;
        int refBaseIndex = refSequence.index(refBaseStartPos);

        for(int i = varReadIndex + 1; i <= min(varReadIndex + indelAltLength, readBases.length - 1); ++i, ++refBaseIndex)
        {
            byte readBase = readBases[i];
            byte refBase = refSequence.Bases[refBaseIndex];

            if(readBase != refBase)
                break;

            if(homology == null)
                homology = new StringBuilder();

            homology.append((char)readBase);
        }

        if(homology == null)
            return null;

        String homologyBases = homology.toString();
        int homologyLength = homologyBases.length();

        if(homologyLength == indelAltLength)
        {
            // continue searching for repeats of the homology to find the point of right-alignment
            boolean matched = true;

            while(matched)
            {
                for(int i = 0; i < indelAltLength; ++i, ++refBaseIndex)
                {
                    byte homologyBase = (byte)homologyBases.charAt(0);
                    byte refBase = refSequence.Bases[refBaseIndex];

                    if(homologyBase != refBase)
                    {
                        matched = false;
                        break;
                    }
                }

                if(matched)
                    homologyLength += indelAltLength;
            }
        }

        return new Microhomology(homologyBases, homologyLength);
    }

    public ReadContext createDelContext(
            final String ref, int refPosition, int readIndex, final byte[] readBases, final RefSequence refSequence)
    {
        /* Routine:
            - test for homology and right-align if required
            - continue in each direction until repeat ends
            - add 2 bases then flank
        */

        // Microhomology microhomology = findHomology(final SimpleVariant variant, final SAMRecord read, int varReadIndex, final RefSequence refSequence)




        /*

        PosCigarBases refPosCigarBases = new PosCigarBases(
            final byte[] bases, final int alignmentStart, final int alignmentEnd,
            final int corePosStart, final int position, final int coreIndexStart,
            final int index, final int coreIndexEnd, final int corePosEnd, final List<CigarElement> cigar)
        */

        return null;
    }
}
