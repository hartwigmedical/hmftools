package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_CHAR_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_CHAR_T;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_TEST_LEN;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.getReadIndexAtReferencePosition;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.common.IndelCoords;

public final class LineUtils
{
    public static int findBaseRepeatCount(final byte[] bases, int index, boolean moveForward, byte testBase)
    {
        // count matching bases until the end of the read
        int matchCount = 0;
        while(index >= 0 && index < bases.length)
        {
            if(bases[index] == testBase)
                ++matchCount;
            else
                break;

            index += moveForward ? 1 : -1;
        }

        return matchCount;
    }

    public static int findLineSequenceCount(final byte[] bases, final int indexStart, final int indexEnd, final byte lineBase)
    {
        if(indexStart < 0 || indexEnd >= bases.length)
            return 0;

        if(indexEnd - indexStart + 1 < LINE_POLY_AT_REQ)
            return 0;

        int lineBaseCount = 0;
        int otherCount = 0;
        for(int i = indexStart; i <= indexEnd; ++i)
        {
            if(bases[i] == lineBase)
            {
                ++lineBaseCount;
            }
            else
            {
                ++otherCount;

                if(otherCount > LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ)
                    break;
            }
        }

        return lineBaseCount >= LINE_POLY_AT_REQ ? lineBaseCount : 0;
    }

    public static int findConsensusLineExtension(final List<Read> reads, final Junction junction)
    {
        // find the poly A/T sequence by first favouring reads heading to the 5' end, ending in non-poly-T bases, taking their median bounds
        Map<Integer,Integer> lengthFrequency = findLineExtensionFrequency(reads, junction, true, true);

        if(!lengthFrequency.isEmpty())
            return calcMedian(lengthFrequency);

        lengthFrequency = findLineExtensionFrequency(reads, junction, false, true);

        if(!lengthFrequency.isEmpty())
            return calcMedian(lengthFrequency);

        lengthFrequency = findLineExtensionFrequency(reads, junction, false, false);

        // take the maximum if from non-terminating sequences
        return lengthFrequency.keySet().stream().mapToInt(x -> x.intValue()).max().orElse(0);
    }

    private static final int MAX_NON_LINE_BASES = LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ;

    @VisibleForTesting
    public static Map<Integer,Integer> findLineExtensionFrequency(
            final List<Read> reads, final Junction junction, boolean fivePrimeOnly, boolean terminatingOnly)
    {
        byte lineBase = junction.isForward() ? LINE_BASE_T : LINE_BASE_A;
        boolean isForward = junction.isForward();

        Map<Integer,Integer> lengthFrequency = Maps.newHashMap();

        for(Read read : reads)
        {
            if(!read.hasLineTail())
                continue;

            if(fivePrimeOnly)
            {
                // only count if moving towards the 5' end of the read
                if(read.positiveStrand() == isForward)
                    continue;
            }

            int readJunctionIndex = getReadIndexAtReferencePosition(read, junction.Position, true);

            if(readJunctionIndex == INVALID_INDEX)
                continue;

            int firstExtensionIndex = readJunctionIndex + (isForward ? 1 : -1);
            int lineExtensionIndex = findLineExtensionEndIndex(read, lineBase, firstExtensionIndex, isForward);
            int lineLength = abs(firstExtensionIndex - lineExtensionIndex) + 1;

            if(lineLength < LINE_POLY_AT_REQ)
                continue;

            if(terminatingOnly)
            {
                // compare LINE length vs soft-clip length
                int softClipLength = isForward ? read.rightClipLength() : read.leftClipLength();

                if(softClipLength == lineLength)
                    continue;
            }

            Integer count = lengthFrequency.get(lineLength);
            lengthFrequency.put(lineLength, count != null ? count + 1 : 1);
        }

        return lengthFrequency;
    }

    public static int findLineExtensionEndIndex(final Read read, final byte lineBase, int extensionIndex, boolean isForward)
    {
        int index = extensionIndex;
        int lineLength = 0;
        int nonLineCount = 0;

        while(index >= 0 && index < read.basesLength())
        {
            if(read.getBases()[index] == lineBase)
            {
                ++lineLength;
            }
            else
            {
                // break on any mismatch once the minimum required LINE site length has been reached
                if(lineLength >= LINE_POLY_AT_REQ)
                    break;

                if(aboveMinQual(read.getBaseQuality()[index]))
                {
                    ++nonLineCount;

                    if(nonLineCount > MAX_NON_LINE_BASES)
                        break;
                }
            }

            index += isForward ? 1 : -1;
        }

        // back up to end of read or last line base
        index += isForward ? -1 : 1;
        return index;
    }

    private static int calcMedian(final Map<Integer,Integer> lengthFrequency)
    {
        if(lengthFrequency.isEmpty())
            return 0;

        int totalEntries = lengthFrequency.values().stream().mapToInt(x -> x.intValue()).sum();

        List<Integer> lengths = Lists.newArrayListWithCapacity(totalEntries);

        for(Map.Entry<Integer,Integer> entry : lengthFrequency.entrySet())
        {
            for(int i = 0; i < entry.getValue(); ++i)
            {
                lengths.add(entry.getKey());
            }
        }

        Collections.sort(lengths);

        int medianIndex = totalEntries / 2;
        return lengths.get(medianIndex);
    }

    public static boolean isLineWithLocalAlignedInsert(final JunctionAssembly assembly)
    {
        // rule out out a line classification if it is supported by reads with inserts with alignment on both sides
        if(assembly.junction().indelBased() || !assembly.hasLineSequence())
            return false;

        int polyAtLength = calcLineSequenceLength(assembly.formJunctionSequence(), assembly.isForwardJunction());
        int maxPostPolyAtLength = 0;

        for(SupportRead read : assembly.support())
        {
            if(read.type() != SupportType.INDEL || read.cachedRead().indelCoords() == null || read.cachedRead().indelCoords().isDelete())
                continue;

            IndelCoords indelCoords = read.cachedRead().indelCoords();

            if(abs(indelCoords.Length - polyAtLength) <= 2)
            {
                int postInsertAlignedLength = 0;

                if(assembly.isForwardJunction())
                    maxPostPolyAtLength = read.alignmentEnd() - indelCoords.PosEnd;
                else
                    maxPostPolyAtLength = indelCoords.PosStart - read.alignmentStart();

                maxPostPolyAtLength = max(maxPostPolyAtLength, postInsertAlignedLength);
            }
        }

        return maxPostPolyAtLength >= MIN_VARIANT_LENGTH;
    }

    public static boolean hasLineSourceSequence(final JunctionAssembly assembly)
    {
        int extIndexStart, extIndexEnd;

        if(assembly.isForwardJunction())
        {
            extIndexStart = assembly.junctionIndex() + 1;
            extIndexEnd = extIndexStart + LINE_POLY_AT_TEST_LEN - 1;
        }
        else
        {
            extIndexEnd = assembly.junctionIndex() - 1;
            extIndexStart = extIndexEnd - LINE_POLY_AT_TEST_LEN + 1;
        }

        byte lineBase = assembly.isForwardJunction() ? LINE_BASE_A : LINE_BASE_T;
        int lineBaseCount = findLineSequenceCount(assembly.bases(), extIndexStart, extIndexEnd, lineBase);
        return lineBaseCount >= LINE_POLY_AT_REQ;
    }

    public static AssemblyLink tryLineSequenceLink(
            final JunctionAssembly first, final JunctionAssembly second, boolean firstReversed, boolean secondReversed)
    {
        if(!first.hasLineSequence() && !second.hasLineSequence())
            return null;

        if(first.hasLineSequence() && second.hasLineSequence()) // also cannot both be LINE insertions - only one will have the poly A/T
            return null;

        // look for a poly A/T match of sufficient length at the ends of each sequence
        String firstExtensionBases = first.formJunctionSequence();
        String firstMatchBases = firstReversed ? Nucleotides.reverseComplementBases(firstExtensionBases) : firstExtensionBases;

        String secondExtensionBases = second.formJunctionSequence();

        if(secondReversed)
            secondExtensionBases = Nucleotides.reverseComplementBases(secondExtensionBases);

        int firstPolyAtLength = calcLineSequenceLength(firstMatchBases, false);
        int secondPolyAtLength = calcLineSequenceLength(secondExtensionBases, true);

        if(firstPolyAtLength < LINE_POLY_AT_REQ || secondPolyAtLength < LINE_POLY_AT_REQ)
            return null;

        String insertedBases = firstExtensionBases;

        int minOverlap = min(firstPolyAtLength, secondPolyAtLength);

        String secondExtraInsertBases = secondExtensionBases.substring(minOverlap);
        insertedBases += secondExtraInsertBases;

        return new AssemblyLink(first, second, LinkType.SPLIT, insertedBases, "");
    }

    public static int calcLineSequenceLength(final String sequence, boolean fromStart)
    {
        int sequenceLength = sequence.length();
        char lineBase = fromStart ? sequence.charAt(0) : sequence.charAt(sequenceLength - 1);

        if(lineBase != LINE_CHAR_A && lineBase != LINE_CHAR_T)
            return 0;

        int count = 1;
        int index = fromStart ? 1 : sequenceLength - 2;
        while(index >= 0 && index < sequenceLength)
        {
            if(sequence.charAt(index) != lineBase)
                break;

            ++count;

            index += fromStart ? 1 : -1;
        }

        return count;
    }
}
