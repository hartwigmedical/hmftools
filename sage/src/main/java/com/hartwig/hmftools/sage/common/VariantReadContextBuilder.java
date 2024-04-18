package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarStringFromElements;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Arrays;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class VariantReadContextBuilder
{
    private final int mFlankSize;

    public VariantReadContextBuilder(int flankSize)
    {
        mFlankSize = flankSize;
    }

    public VariantReadContext createSnvMnvContext(
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
                contextReadBases, readCigarInfo.Cigar, coreIndexStart, readVarIndex, coreIndexEnd,
                null, repeatBoundaryInfo.MaxRepeat);
    }

    public VariantReadContext createIndelContext(
            final SimpleVariant variant, final SAMRecord read, int varReadIndex, final RefSequence refSequence)
    {
        /* Routine:
            - test for homology and right-align if required
            - continue in each direction until repeat ends
            - add 2 bases then flank
        */

        // note that for indels only the 2 ref bases on each side are part of the core
        int readCoreStart = varReadIndex - MIN_CORE_DISTANCE + 1;

        int readCoreEnd = varReadIndex + (variant.isInsert() ? variant.indelLength() + 1 : 1) + MIN_CORE_DISTANCE - 1;

        final byte[] readBases = read.getReadBases();

        Microhomology homology = findHomology(variant, read, varReadIndex, refSequence);

        if(homology != null)
            readCoreEnd += homology.Length;

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
                contextReadBases, readCigarInfo.Cigar, coreIndexStart, readVarIndex, coreIndexEnd, homology, repeatBoundaryInfo.MaxRepeat);
    }

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

        public String toString()
        {
            return format("core(%d-%d) repeat(%s)", CoreStart, CoreEnd, MaxRepeat != null ? MaxRepeat : "");
        }
    }

    // set an initial search length long enough to find a min count of the longest repeat
    private static final int REPEAT_SEARCH_LENGTH = MAX_REPEAT_LENGTH * MIN_REPEAT_COUNT;

    private RepeatBoundaryInfo findRepeatBoundaries(int lowerRefIndex, int upperRefIndex, final byte[] readBases)
    {
        int repeatSearchStart = max(lowerRefIndex - REPEAT_SEARCH_LENGTH, 0);

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

        // stored for the scenario where an indel in the flanks pushes out the alignment beyond the standard flank length
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

        public String toString()
        {
            return format("%s align(%d-%d) flankIndex(%d-%d)",
                    cigarStringFromElements(Cigar), AlignmentStart, AlignmentEnd, FlankIndexStart, FlankIndexEnd);
        }
    }

    private ReadCigarInfo buildReadCigar(final SAMRecord read, int indexStart, int indexEnd)
    {
        // builds a cigar around the specified index boundaries and calculates the corresponding alignment positions
        List<CigarElement> cigar = Lists.newArrayList();

        int readIndex = 0;
        int refPosition = read.getAlignmentStart();
        int alignmentStart = 0;
        int alignmentEnd = 0;

        int finalIndexStart = indexStart;
        int finalIndexEnd = indexEnd;

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
                    if(element.getOperator() == I)
                    {
                        // handles an insert that pushes the alignment out - always take the prior alignment base and reduce index start
                        // eg looking to find the alignment boundary for index start 12 for 10M5I1350M, alignment start == 100
                        // so at the insert element, read index = 10, ref pos = 110 (pointing at next ref base)
                        int extraIndexStart = indexStart - readIndex + 1;
                        cigar.add(0, new CigarElement(1, M));

                        finalIndexStart -= extraIndexStart;
                        alignmentStart = max(refPosition - 1, read.getAlignmentStart());
                    }
                    else if(element.getOperator() == D)
                    {
                        alignmentStart = max(refPosition - 1, read.getAlignmentStart());
                    }
                    else
                    {
                        alignmentStart = refPosition + (indexStart - readIndex);
                    }
                }

                if(alignmentEnd == 0 && elementEndIndex >= indexEnd)
                {
                    if(element.getOperator() == I)
                    {
                        // similar extension to the above
                        int extraIndexEnd = elementEndIndex + 1 - indexEnd;
                        cigar.add(new CigarElement(1, M));

                        finalIndexEnd += extraIndexEnd;

                        alignmentEnd = refPosition; // already pointing at the next M (aligned) base
                    }
                    else if(element.getOperator() == D)
                    {
                        alignmentEnd = refPosition + element.getLength();
                    }
                    else
                    {
                        alignmentEnd = min(refPosition + (indexEnd - readIndex), read.getAlignmentEnd());
                    }
                }
            }

            if(element.getOperator().consumesReadBases())
                readIndex += element.getLength();

            if(element.getOperator().consumesReferenceBases())
                refPosition += element.getLength();

            if(readIndex > indexEnd)
                break;
        }

        return new ReadCigarInfo(cigar, alignmentStart, alignmentEnd, finalIndexStart, finalIndexEnd);
    }

    public static Microhomology findHomology(final SimpleVariant variant, final SAMRecord read, int varReadIndex, final RefSequence refSequence)
    {
        if(!variant.isIndel())
            return null;

        int indelAltLength = abs(variant.indelLength());

        StringBuilder homology = null;
        String indelBases = variant.isInsert() ? variant.alt().substring(1) : variant.ref().substring(1);
        // final byte[] readBases = read.getReadBases();
        int refBaseStartPos = variant.isInsert() ? variant.Position + 1 : variant.Position + indelAltLength + 1;
        int refBaseIndex = refSequence.index(refBaseStartPos);

        for(int i = 0; i < indelBases.length(); ++i, ++refBaseIndex)
        {
            byte readBase = (byte)indelBases.charAt(i);
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
}
