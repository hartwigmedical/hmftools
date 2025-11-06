package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.sage.SageConfig.SEQUENCING_TYPE;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConfig.isIllumina;
import static com.hartwig.hmftools.sage.SageConfig.isUltima;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.common.SageVariant.isLongInsert;
import static com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender.extendUltimaCore;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.ultimaLongRepeatFilter;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Collections;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.seqtech.IlluminaArtefactContext;
import com.hartwig.hmftools.sage.seqtech.UltimaCoreExtender;

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
            int minFlankLength = isUltima() && SageConfig.AppendMode ? 1 : mFlankSize;
            if(readContext.leftFlankLength() < minFlankLength || readContext.rightFlankLength() < minFlankLength)
                return null;

            if(!aboveMinBaseQual(readContext, read, varReadIndex))
                return null;

            // enforce ref base padding for trinucleotide generation
            if(readContext.variantRefIndex() <= 0 || readContext.variantRefIndex() >= readContext.refBases().length() - 1)
                return null;

            // ref bases are extended for long inserts for partial ref matching
            if(isLongInsert(variant))
            {
                byte[] extendedRefBases = refSequence.baseRange(
                        readContext.CorePositionStart - readContext.leftFlankLength(),
                        readContext.CorePositionEnd + readContext.rightFlankLength() + 1);
                readContext.setExtendedRefBases(extendedRefBases);
            }

            // set max ref repeat for use in MSI calcs and VCF output
            setMaxRefRepeat(readContext);

            if(isIllumina())
                readContext.setArtefactContext(IlluminaArtefactContext.buildContext(readContext));

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

    @VisibleForTesting
    public VariantReadContext buildContext(
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
        Microhomology refHomology = null;
        SoftClipReadAdjustment softClipReadAdjustment = null;

        if(variant.isIndel())
        {
            readCoreStart = varIndexInRead - MIN_CORE_DISTANCE + 1;
            readCoreEnd = varIndexInRead + (variant.isInsert() ? variant.indelLength() + 1 : 1) + MIN_CORE_DISTANCE - 1;
            int refBaseLength = read.getReadBases().length - varIndexInRead + (variant.isDelete() ? variant.indelLengthAbs() : 0);
            byte[] homologyRefBases = refSequence.baseRange(variant.position(), variant.position() + refBaseLength);

            softClipReadAdjustment = checkIndelSoftClipAdjustment(read, variant, varIndexInRead);

            if(softClipReadAdjustment != null)
            {
                homology = Microhomology.findHomology(variant, homologyRefBases, 0, false);
            }
            else
            {
                homology = Microhomology.findHomology(variant, read.getReadBases(), varIndexInRead, variant.isInsert());
                refHomology = Microhomology.findHomology(variant, homologyRefBases, 0, variant.isDelete());
            }

            if(refHomology != null && (homology == null || refHomology.Length  > homology.Length))
                homology = refHomology;
            if(homology != null)
                readCoreEnd += homology.Length;

            if(readCoreEnd + mFlankSize >= read.getReadBases().length && homology != null && SageConfig.AppendMode && isLongInsert(variant))
            {
                // most likely a long insert aligned from softclip - assume the old core is correct
                readCoreEnd = read.getReadBases().length - mFlankSize - 1;
            }
        }
        else
        {
            readCoreStart = varIndexInRead - MIN_CORE_DISTANCE;
            readCoreEnd = varIndexInRead + variant.Alt.length() - 1 + MIN_CORE_DISTANCE;
        }

        if(readCoreStart < 0 || readCoreEnd >= read.getReadBases().length)
            return null;

        if(isUltima() && ultimaLongRepeatFilter(variant, read, varIndexInRead, homology))
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

        // for ultima we expand core so that homopolymers are not cut off in the read or the ref
        if(isUltima())
        {
            // long inserts with S elements won't work properly in append mode with core extension
            boolean skippableLongInsert = SageConfig.AppendMode && isLongInsert(variant);

            UltimaCoreExtender.UltimaCoreInfo ultimaCoreInfo = extendUltimaCore(
                    read.getReadBases(), refSequence,
                    softClipReadAdjustment != null ? softClipReadAdjustment.AlignmentStart : read.getAlignmentStart(),
                    softClipReadAdjustment != null ? softClipReadAdjustment.ConvertedCigar : read.getCigar().getCigarElements(),
                    readCigarInfo, mFlankSize, SageConfig.AppendMode);

            if(!skippableLongInsert && (ultimaCoreInfo == null || !ultimaCoreInfo.CigarInfo.isValid()))
                return null;

            if(ultimaCoreInfo != null)
            {
                readCoreStart = ultimaCoreInfo.ReadCoreStart;
                readCoreEnd = ultimaCoreInfo.ReadCoreEnd;
                readCigarInfo = ultimaCoreInfo.CigarInfo;
            }
        }

        int readPositionStart = readCigarInfo.ReadAlignmentStart; // may have been adjusted
        readFlankStart = readCigarInfo.FlankIndexStart;
        readFlankEnd = readCigarInfo.FlankIndexEnd;

        if(readFlankStart < 0 || readFlankEnd >= readBases.length)
            return null;

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

    private static boolean aboveMinBaseQual(final VariantReadContext readContext, final SAMRecord read, final int varReadIndex)
    {
        int indexStart = varReadIndex - readContext.leftCoreLength() - readContext.leftFlankLength();
        int indexEnd = varReadIndex + readContext.rightCoreLength() + readContext.rightFlankLength();

        if(indexStart < 0 || indexEnd >= read.getBaseQualities().length)
            return false;

        for(int i = indexStart; i <= indexEnd; ++i)
        {
            if(BaseQualAdjustment.isUncertainBaseQual(read.getBaseQualities()[i])
            || BaseQualAdjustment.isMediumBaseQual(read.getBaseQualities()[i], SEQUENCING_TYPE))
                return false;
        }

        return true;
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

    public static void setMaxRefRepeat(final VariantReadContext readContext)
    {
        // set max ref repeat for use in MSI calcs and VCF output
        int refRepeatIndex = readContext.variant().isIndel() ? readContext.variantRefIndex() + 1 : readContext.variantRefIndex();

        RepeatInfo maxRefRepeat = RepeatInfo.findMaxRepeat(
                readContext.RefBases, refRepeatIndex, refRepeatIndex, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT + 1,
                false, refRepeatIndex);

        readContext.setRefMaxRepeat(maxRefRepeat);
    }
}
