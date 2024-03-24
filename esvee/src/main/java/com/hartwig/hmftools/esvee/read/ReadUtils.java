package com.hartwig.hmftools.esvee.read;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.common.SvConstants.DISCORDANT_FRAGMENT_LENGTH;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.AssemblyConstants;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public final class ReadUtils
{
    public static String readToString(final SAMRecord read)
    {
        return format("id(%s) coords(%s:%d-%d) cigar(%s) mate(%s:%d) flags(%d)",
                read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(),
                read.getCigarString(), read.getMateReferenceName(), read.getMateAlignmentStart(), read.getFlags());
    }

    public static boolean isDiscordant(final Read read)
    {
        return isDiscordant(read, DISCORDANT_FRAGMENT_LENGTH);
    }

    public static boolean isDiscordant(final Read read, final int discordantPairFragmentLength)
    {
        // FIXME: share method from SvUtils and/or SvPrep
        if(read.isMateUnmapped())
            return false;

        if(!read.chromosome().equals(read.mateChromosome())) // not strictly correct for supplementaries, since needs to check the primary
            return true;

        if(read.positiveStrand() == read.matePositiveStrand())
            return true;

        int fragmentSize = abs(read.insertSize());

        return fragmentSize == 0 || fragmentSize >= discordantPairFragmentLength;
    }

    public static final int INVALID_INDEX = -1;

    public static int getReadIndexAtReferencePosition(final Read read, final int refPosition, boolean allowExtrapolation)
    {
        // finds the read index given a reference position, and extrapolates outwards from alignments as required

        // for indel reads, use the implied alignment and unclipped positions from each direction
        int alignmentStart = allowExtrapolation && read.indelImpliedAlignmentStart() > 0 ?
                read.indelImpliedAlignmentStart() : read.alignmentStart();

        int alignmentEnd = allowExtrapolation && read.indelImpliedAlignmentEnd() > 0 ?
                read.indelImpliedAlignmentEnd() : read.alignmentEnd();

        if(refPosition <= alignmentStart)
        {
            if(!allowExtrapolation && refPosition < alignmentStart)
                return INVALID_INDEX;

            int baseDiff = alignmentStart - refPosition;
            int softClipBases = alignmentStart - read.unclippedStart();
            return baseDiff <= softClipBases ? softClipBases - baseDiff : INVALID_INDEX;
        }
        else if(refPosition >= alignmentEnd)
        {
            if(!allowExtrapolation && refPosition > alignmentEnd)
                return INVALID_INDEX;

            int baseDiff = refPosition - alignmentEnd;
            int softClipBases = read.unclippedEnd() - alignmentEnd;
            return baseDiff <= softClipBases ? read.basesLength() - (softClipBases - baseDiff) - 1 : INVALID_INDEX;
        }

        // cannot use standard method since CIGAR and coords may have been adjusted
        int readIndex = 0;
        int currentPos = read.alignmentStart();
        for(CigarElement element : read.cigarElements())
        {
            if(!element.getOperator().consumesReferenceBases())
            {
                readIndex += element.getLength();
                continue;
            }

            if(currentPos == refPosition)
                break;

            if(!element.getOperator().consumesReadBases())
            {
                // for a D or N where the position is inside it, return the read index for the start of the element
                if(refPosition >= currentPos && refPosition < currentPos + element.getLength())
                    return readIndex - 1;

                currentPos += element.getLength();
            }
            else
            {
                // pos = 100, element = 10M, covering pos 100-109, read index 4 (say after 4S), ref pos at last base of element = 109
                if(refPosition >= currentPos && refPosition < currentPos + element.getLength())
                    return readIndex + refPosition - currentPos;

                currentPos += element.getLength();
                readIndex += element.getLength();
            }
        }

        return readIndex;
    }

    public static int avgBaseQuality(final Read read) { return avgBaseQuality(read.getBaseQuality(), 0, read.basesLength() - 1); }

    public static int avgBaseQuality(final byte[] baseQualities, final int startIndex, final int endIndex)
    {
        if(startIndex > endIndex || startIndex < 0 || endIndex >= baseQualities.length)
            return -1;

        int qualitySum = 0;
        for(int i = startIndex; i <= endIndex; i++)
        {
            qualitySum += baseQualities[i];
        }

        return qualitySum / (endIndex - startIndex + 1);
    }

    public static void copyArray(final byte[] source, final byte[] dest, final int sourceIndexStart, final int destIndexStart)
    {
        int d = destIndexStart;
        for(int s = sourceIndexStart; s < source.length && d < dest.length; ++s, ++d)
        {
            dest[d] = source[s];
        }
    }

    public static byte[] copyArray(final byte[] source)
    {
        byte[] dest = new byte[source.length];

        for(int i = 0; i < source.length; ++i)
        {
            dest[i] = source[i];
        }

        return dest;
    }

    public static byte[] subsetArray(final byte[] source, final int startIndex, final int endIndex)
    {
        byte[] dest = new byte[endIndex - startIndex + 1];

        int newIndex = 0;
        for(int index = startIndex; index <= endIndex; ++index, ++newIndex)
        {
            dest[newIndex] = source[index];
        }

        return dest;
    }

    public static int[] copyArray(final int[] source)
    {
        int[] dest = new int[source.length];

        for(int i = 0; i < source.length; ++i)
        {
            dest[i] = source[i];
        }

        return dest;
    }

    public static byte[] addByteArray(final byte[] first, final byte[] second)
    {
        byte[] combined = new byte[first.length + second.length];

        for(int i = 0; i < first.length; ++i)
        {
            combined[i] = first[i];
        }

        for(int i = 0; i < second.length; ++i)
        {
            combined[first.length + i] = second[i];
        }

        return combined;
    }

    public static byte[] reverseBytes(final byte[] bases)
    {
        String reversed = Nucleotides.reverseComplementBases(new String(bases));
        return reversed.getBytes();
    }

    public static void initialise(final byte[] array, final byte value)
    {
        for(int i = 0; i < array.length; ++i)
        {
            array[i] = value;
        }
    }

    public static void initialise(final int[] array, final int value)
    {
        for(int i = 0; i < array.length; ++i)
        {
            array[i] = value;
        }
    }
}
