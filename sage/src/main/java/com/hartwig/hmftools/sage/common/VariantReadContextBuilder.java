package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.common.SimpleVariant.isLongInsert;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.sage.quality.ArtefactContext;

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

            // ref bases are extended for long inserts for partial ref matching
            if(isLongInsert(variant))
            {
                readContext.setExtendedRefBases(refSequence.positionBases(
                        readContext.CorePositionEnd + 1, readContext.CorePositionEnd + 11));
            }

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
        SoftClipReadAdjustment softClipReadAdjustment = null;

        if(variant.isIndel())
        {
            readCoreStart = varIndexInRead - MIN_CORE_DISTANCE + 1;
            readCoreEnd = varIndexInRead + (variant.isInsert() ? variant.indelLength() + 1 : 1) + MIN_CORE_DISTANCE - 1;

            softClipReadAdjustment = checkIndelSoftClipAdjustment(read, variant, varIndexInRead);

            if(softClipReadAdjustment != null)
            {
                int refBaseLength = read.getReadBases().length - varIndexInRead;
                byte[] homologyRefBases = refSequence.baseRange(variant.position(), variant.position() + refBaseLength);
                homology = Microhomology.findHomology(variant, homologyRefBases, 0, false);
            }
            else
            {
                homology = Microhomology.findHomology(variant, read.getReadBases(), varIndexInRead, true);
            }

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

        readCoreStart = min(readCoreStart, repeatBoundaries.LowerIndex);
        readCoreEnd = max(readCoreEnd, repeatBoundaries.UpperIndex);

        int readFlankStart = readCoreStart - mFlankSize;
        int readFlankEnd = readCoreEnd + mFlankSize;

        if(readFlankStart < 0 || readFlankEnd >= read.getReadBases().length)
            return null;

        // build a CIGAR from the read to cover this range
        ReadCigarInfo readCigarInfo = ReadCigarInfo.buildReadCigar(
                softClipReadAdjustment != null ? softClipReadAdjustment.AlignmentStart : read.getAlignmentStart(),
                softClipReadAdjustment != null ? softClipReadAdjustment.ConvertedCigar : read.getCigar().getCigarElements(),
                readFlankStart, readCoreStart, readCoreEnd, readFlankEnd);

        if(readCigarInfo == null || !readCigarInfo.isValid())
            return null;

        int readPositionStart = readCigarInfo.ReadAlignmentStart; // may have been adjusted
        readFlankStart = readCigarInfo.FlankIndexStart;
        readFlankEnd = readCigarInfo.FlankIndexEnd;

        int corePositionStart = readCigarInfo.CorePositionStart;
        int corePositionEnd = readCigarInfo.CorePositionEnd;

        int alignmentStart = max(readPositionStart, readCigarInfo.FlankPositionStart);
        int alignmentEnd = min(read.getAlignmentEnd(), readCigarInfo.FlankPositionEnd);

        byte[] contextReadBases = Arrays.subsetArray(readBases, readFlankStart, readFlankEnd);

        int readContextOffset = readFlankStart;
        int coreIndexStart = readCoreStart - readContextOffset;
        int readVarIndex = varIndexInRead - readContextOffset;
        int coreIndexEnd = readCoreEnd - readContextOffset;

        // ref bases are the core width around the variant's position
        byte[] refBases = refSequence.baseRange(corePositionStart, corePositionEnd);

        RepeatInfo maxRepeat = null;
        List<RepeatInfo> allRepeats;

        if(repeatBoundaries.MaxRepeat != null)
        {
            maxRepeat = new RepeatInfo(
                    repeatBoundaries.MaxRepeat.Index - readContextOffset,
                    repeatBoundaries.MaxRepeat.Bases, repeatBoundaries.MaxRepeat.Count);
        }

        if(repeatBoundaries.AllRepeats != null)
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
                variant, alignmentStart, alignmentEnd, refBases, contextReadBases, readCigarInfo.Cigar, coreIndexStart,
                readVarIndex, coreIndexEnd, homology, maxRepeat, allRepeats, corePositionStart, corePositionEnd);
    }

    private class SoftClipReadAdjustment
    {
        public final boolean IsLeft;
        public final int AlignmentStart;
        public final List<CigarElement> ConvertedCigar;

        public SoftClipReadAdjustment(final boolean isLeft, final int alignmentStart, final List<CigarElement> convertedCigar)
        {
            IsLeft = isLeft;
            AlignmentStart = alignmentStart;
            ConvertedCigar = convertedCigar;
        }
    }

    private SoftClipReadAdjustment checkIndelSoftClipAdjustment(final SAMRecord read, final SimpleVariant variant, int varReadIndex)
    {
        if(!isLongInsert(variant))
            return null;

        if(read.getCigar().getCigarElements().get(0).getOperator() == S
        && varReadIndex < read.getCigar().getCigarElements().get(0).getLength())
        {
            // correct inserts in soft-clips before building the cigar
            List<CigarElement> origReadCigar = read.getCigar().getCigarElements();

            // turn soft-clips into inserts when they originate in a soft-clip
            int leftSoftClipLength = origReadCigar.get(0).getLength();

            List<CigarElement> convertedCigar = Lists.newArrayList(origReadCigar);
            convertedCigar.remove(0);

            convertedCigar.add(0, new CigarElement(varReadIndex + 1, M));
            convertedCigar.add(1, new CigarElement(variant.indelLength(), I));

            int remainingAlignedBases = leftSoftClipLength - (varReadIndex + 1) - variant.indelLength();

            if(remainingAlignedBases > 0)
            {
                int postInsertAlignedBases = convertedCigar.get(2).getLength();
                convertedCigar.remove(2);
                convertedCigar.add(2, new CigarElement(remainingAlignedBases + postInsertAlignedBases, M));
            }

            int convertedStartPosition = variant.position() - varReadIndex;

            return new SoftClipReadAdjustment(true, convertedStartPosition, convertedCigar);
        }

        int lastCigarIndex = read.getCigar().getCigarElements().size() - 1;
        if(read.getCigar().getCigarElements().get(lastCigarIndex).getOperator() == S)
        {
            List<CigarElement> origReadCigar = read.getCigar().getCigarElements();
            int rightSoftClipLength = origReadCigar.get(lastCigarIndex).getLength();
            int rightSoftClipStartIndex = read.getReadBases().length - rightSoftClipLength;

            if(varReadIndex >= rightSoftClipStartIndex)
            {
                List<CigarElement> convertedCigar = Lists.newArrayList(origReadCigar);
                convertedCigar.remove(lastCigarIndex);

                convertedCigar.add(new CigarElement(variant.indelLength(), I));

                int remainingAlignedBases = rightSoftClipLength - variant.indelLength();
                convertedCigar.add(new CigarElement(remainingAlignedBases, M));

                return new SoftClipReadAdjustment(false, read.getAlignmentStart(), convertedCigar);
            }
        }

        return null;
    }

    private RepeatBoundaries findRepeatBoundaries(int readCoreStart, int readCoreEnd, final byte[] readBases)
    {
        return RepeatBoundaries.findRepeatBoundaries(readBases, readCoreStart, readCoreEnd, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);
    }
}
