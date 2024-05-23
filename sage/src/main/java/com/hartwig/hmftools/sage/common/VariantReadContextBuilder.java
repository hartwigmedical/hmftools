package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;

import java.util.Collections;
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

            homology = Microhomology.findHomology(variant, read, varIndexInRead);

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

        if(repeatBoundaries != null)
        {
            readCoreStart = min(readCoreStart, repeatBoundaries.LowerIndex);
            readCoreEnd = max(readCoreEnd, repeatBoundaries.UpperIndex);
        }

        int readFlankStart = readCoreStart - mFlankSize;
        int readFlankEnd = readCoreEnd + mFlankSize;

        if(readFlankStart < 0 || readFlankEnd >= read.getReadBases().length)
            return null;

        // build a CIGAR from the read to cover this range
        ReadCigarInfo readCigarInfo = ReadCigarInfo.buildReadCigar(read, readFlankStart, readCoreStart, readCoreEnd, readFlankEnd);

        if(readCigarInfo == null)
            return null;

        readFlankStart = readCigarInfo.FlankIndexStart;
        readFlankEnd = readCigarInfo.FlankIndexEnd;

        int corePositionStart = readCigarInfo.CorePositionStart;
        int corePositionEnd = readCigarInfo.CorePositionEnd;

        int alignmentStart = max(read.getAlignmentStart(), readCigarInfo.FlankPositionStart);
        int alignmentEnd = min(read.getAlignmentEnd(), readCigarInfo.FlankPositionEnd);

        byte[] contextReadBases = Arrays.subsetArray(readBases, readFlankStart, readFlankEnd);

        int readContextOffset = readFlankStart;
        int coreIndexStart = readCoreStart - readContextOffset;
        int readVarIndex = varIndexInRead - readContextOffset;
        int coreIndexEnd = readCoreEnd - readContextOffset;

        // ref bases are the core width around the variant's position
        byte[] refBases = refSequence.baseRange(corePositionStart, corePositionEnd);

        int altIndexLower = determineLowerAltIndex(variant, readVarIndex);
        int refIndex = readVarIndex - coreIndexStart;
        int altIndexUpper = determineUpperAltIndex(variant, contextReadBases, refBases, readVarIndex, refIndex, coreIndexEnd);

        RepeatInfo maxRepeat = repeatBoundaries != null ? repeatBoundaries.MaxRepeat : null;

        List<RepeatInfo> allRepeats;

        if(repeatBoundaries != null)
        {
            allRepeats = Lists.newArrayListWithCapacity(repeatBoundaries.AllRepeats.size());
            int readFlankOffset = readFlankStart;
            repeatBoundaries.AllRepeats.forEach(x -> allRepeats.add(new RepeatInfo(
                    x.Index - readFlankOffset, x.Bases, x.Count)));
        }
        else
        {
            allRepeats = Collections.emptyList();
        }

        return new VariantReadContext(
                variant, alignmentStart, alignmentEnd, refBases, contextReadBases, readCigarInfo.Cigar, coreIndexStart, readVarIndex,
                coreIndexEnd, homology, maxRepeat, allRepeats, altIndexLower, altIndexUpper, corePositionStart, corePositionEnd);
    }

    private RepeatBoundaries findRepeatBoundaries(int readCoreStart, int readCoreEnd, final byte[] readBases)
    {
        return RepeatBoundaries.findRepeatBoundaries(readBases, readCoreStart, readCoreEnd, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);
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

    public static int determineLowerAltIndex(final SimpleVariant variant, final int varReadIndex)
    {
        if(!variant.isInsert())
            return varReadIndex;

        // return the last inserted base for insert, since this will be the first point of difference on the lower side
        return varReadIndex + variant.indelLength();
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
