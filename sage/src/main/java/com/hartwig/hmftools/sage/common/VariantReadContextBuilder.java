package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
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
import com.hartwig.hmftools.sage.quality.ArtefactContext;

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

            homology = Microhomology.findHomology(variant, read.getReadBases(), varIndexInRead);

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

            /*
            // expand alt boundaries to cover any repeats on that side
            int minRepeatIndex = allRepeats.get(0).Index - 1;
            int maxRepeatIndex = allRepeats.get(0).endIndex() + 1;

            for(int i = 1; i < allRepeats.size(); ++i)
            {
                minRepeatIndex = min(minRepeatIndex, allRepeats.get(i).Index - 1);
                maxRepeatIndex = max(maxRepeatIndex, allRepeats.get(i).endIndex() + 1);
            }

            altIndexLower = min(altIndexLower, minRepeatIndex);
            altIndexUpper = max(altIndexUpper, maxRepeatIndex);
            */
        }
        else
        {
            allRepeats = Collections.emptyList();
        }

        return new VariantReadContext(
                variant, alignmentStart, alignmentEnd, refBases, contextReadBases, readCigarInfo.Cigar, coreIndexStart, readVarIndex,
                coreIndexEnd, homology, maxRepeat, allRepeats, corePositionStart, corePositionEnd);
    }

    private RepeatBoundaries findRepeatBoundaries(int readCoreStart, int readCoreEnd, final byte[] readBases)
    {
        return RepeatBoundaries.findRepeatBoundaries(readBases, readCoreStart, readCoreEnd, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);
    }
}
