package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.SvConstants.BAM_HEADER_SAMPLE_ID_TAG;

import com.hartwig.hmftools.esvee.SvConstants;

public final class ReadUtils
{
    public static boolean isDiscordant(final Read read)
    {
        return isDiscordant(read, SvConstants.DISCORDANT_FRAGMENT_LENGTH);
    }

    public static boolean isDiscordant(final Read read, final int discordantPairFragmentLength)
    {
        if(!read.isMateMapped())
            return false;

        if(!read.chromosome().equals(read.mateChromosome()))
            return true;

        if(read.positiveStrand() == read.matePositiveStrand())
            return true;

        int fragmentSize = read.insertSize();

        return fragmentSize == 0 || fragmentSize >= discordantPairFragmentLength;
    }

    public static int getAvgBaseQuality(final Read read, final int startPosition, final int length)
    {
        final byte[] baseQualities = read.getBaseQuality();
        final int startIndex = startPosition - 1;
        final int endIndex = Math.min(startIndex + length, baseQualities.length);

        int qualitySum = 0;
        for(int i = startIndex; i < endIndex; i++)
            qualitySum += baseQualities[i];
        return qualitySum / length;
    }

    public static int getReadPositionAtReferencePosition(final Read read, final int position)
    {
        // CHECK: is this the same implementation as below? obviously needs to adjust for any changes to cigar
        return read.bamRecord().getReadPositionAtReferencePosition(position);
    }

    public static int getReadIndexAtReferencePosition(final Read read, final int position)
    {
        // finds the read index given a reference position, and extrapolates outwards from alignments as required

        return read.bamRecord().getReadPositionAtReferencePosition(position);
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
