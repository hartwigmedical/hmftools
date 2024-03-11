package com.hartwig.hmftools.esvee.read;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.SvConstants.BAM_HEADER_SAMPLE_ID_TAG;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.SvConstants;

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
        return isDiscordant(read, SvConstants.DISCORDANT_FRAGMENT_LENGTH);
    }

    public static boolean isDiscordant(final Read read, final int discordantPairFragmentLength)
    {
        // FIXME: share method from SvUtils and/or SvPrep
        if(read.isMateUnmapped())
            return false;

        if(!read.chromosome().equals(read.mateChromosome()))
            return true;

        if(read.positiveStrand() == read.matePositiveStrand())
            return true;

        int fragmentSize = abs(read.insertSize());

        return fragmentSize == 0 || fragmentSize >= discordantPairFragmentLength;
    }

    public static double avgBaseQuality(final Read read)
    {
        byte[] baseQualities = read.getBaseQuality();

        int qualitySum = 0;
        for(int i = 0; i < baseQualities.length; i++)
        {
            qualitySum += baseQualities[i];
        }

        return qualitySum / (double)baseQualities.length;
    }

    public static int avgBaseQuality(final Read read, final int startPosition, final int length)
    {
        byte[] baseQualities = read.getBaseQuality();
        int startIndex = startPosition - 1;
        int endIndex = Math.min(startIndex + length, baseQualities.length);

        int qualitySum = 0;
        for(int i = startIndex; i < endIndex; i++)
        {
            qualitySum += baseQualities[i];
        }

        return qualitySum / length;
    }

    @Deprecated
    public static int getReadPositionAtReferencePosition(final Read read, final int position)
    {
        // CHECK: is this the same implementation as below? obviously needs to adjust for any changes to cigar
        return read.bamRecord().getReadPositionAtReferencePosition(position);
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

    /*
    public static int getReadPositionAtReferencePosition(final Read read, final int position)
    {
        if(position <= 0)
            return 0;

        for(Alignment alignmentBlock : read.getAlignmentBlocks())
        {
            final int end = alignmentBlock.ReferenceStartPosition + alignmentBlock.Length - 1;
            if(end >= position)
            {
                if(position < alignmentBlock.ReferenceStartPosition)
                    return 0;
                else
                    return position - alignmentBlock.ReferenceStartPosition + alignmentBlock.SequenceStartPosition;
            }
        }
        return 0;
    }
    */
    
    public static Read flipRead(final Read read)
    {
        // TO-DO - either skip doing this an just register a reversed read with assembly / support, or mark the read or copy it
        return read;
    }

    public static boolean isGermline(final Read read, final String refSampleId)
    {
        return read.bamRecord().getHeader().getAttribute(BAM_HEADER_SAMPLE_ID_TAG).equals(refSampleId);
    }

    /*
    default MutableRecord trimLeft(final int count)
    {
        final MutableRecord clone = copyRecord();
        clone.setBases(Arrays.copyOfRange(getBases(), count, getLength()),
                Arrays.copyOfRange(getBaseQuality(), count, getLength()));
        if(!clone.isUnmapped())
        {
            clone.setCigar(CigarUtils.trimLeft(clone.getCigar(), count));
            final int alignmentMove = Math.max(0, count - leftSoftClipLength(getCigar()));
            clone.setAlignmentStart(getAlignmentStart() + alignmentMove);
        }

        return clone;
    }

    default MutableRecord trimRight(final int count)
    {
        final MutableRecord clone = copyRecord();
        clone.setBases(Arrays.copyOfRange(getBases(), 0, getLength() - count),
                Arrays.copyOfRange(getBaseQuality(), 0, getLength() - count));
        if(!clone.isUnmapped())
            clone.setCigar(CigarUtils.trimRight(clone.getCigar(), count));

        return clone;
    }

    default MutableRecord flipRecord()
    {
        final byte[] readBases = getBases().clone();
        SequenceUtil.reverseComplement(readBases);
        final byte[] newQuals = getBaseQuality().clone();
        SequenceUtil.reverseQualities(newQuals);

        final MutableRecord flipped = copyRecord();
        flipped.setBases(readBases, newQuals);
        flipped.setPositiveStrand(!isPositiveStrand());

        final var cigarElements = new ArrayList<>(flipped.getCigar().getCigarElements());
        Collections.reverse(cigarElements);
        flipped.setCigar(new Cigar(cigarElements));

        return flipped;
    }
    */

}
