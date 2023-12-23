package com.hartwig.hmftools.esvee.models;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import com.hartwig.hmftools.esvee.sam.CigarUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.SequenceUtil;

public interface MutableRecord extends IRecord
{
    void setBases(final byte[] bases, final byte[] quals);

    void setUnmapped();

    void setChromosome(final String chromosome);

    void setAlignmentStart(final int position);

    void setMappingQuality(final int mapQ);

    void setCigar(final Cigar cigar);

    void setCigar(final String cigar);

    void setMateChromosome(final String chromosome);

    /** 0 for no alignment */
    void setMateAlignmentStart(final int position);

    void setMateUnmapped();

    void setPositiveStrand(final boolean isPositiveStrand);

    void setMatePositiveStrand(final boolean isPositiveStrand);

    <T> void setAttribute(final String name, final T value);

    void setIsGermline(final boolean isGermline);

    void setReadPairedFlag(final boolean isPaired);

    void setIsFirstOfPair(final boolean isFirstOfPair);

    MutableRecord copy();

    default MutableRecord trimLeft(final int count)
    {
        final MutableRecord clone = copy();
        clone.setBases(Arrays.copyOfRange(getBases(), count, getLength()),
                Arrays.copyOfRange(getBaseQuality(), count, getLength()));
        if (!clone.isUnmapped())
        {
            clone.setCigar(CigarUtils.trimLeft(clone.getCigar(), count));
            final int alignmentMove = Math.max(0, count - CigarUtils.unmappedStartBases(getCigar()));
            clone.setAlignmentStart(getAlignmentStart() + alignmentMove);
        }

        return clone;
    }

    default MutableRecord trimRight(final int count)
    {
        final MutableRecord clone = copy();
        clone.setBases(Arrays.copyOfRange(getBases(), 0, getLength() - count),
                Arrays.copyOfRange(getBaseQuality(), 0, getLength() - count));
        if (!clone.isUnmapped())
            clone.setCigar(CigarUtils.trimRight(clone.getCigar(), count));

        return clone;
    }

    default MutableRecord flipStrand()
    {
        final byte[] readBases = getBases().clone();
        SequenceUtil.reverseComplement(readBases);
        final byte[] newQuals = getBaseQuality().clone();
        SequenceUtil.reverseQualities(newQuals);

        final MutableRecord flipped = copy();
        flipped.setBases(readBases, newQuals);
        flipped.setPositiveStrand(!isPositiveStrand());

        final var cigarElements = new ArrayList<>(flipped.getCigar().getCigarElements());
        Collections.reverse(cigarElements);
        flipped.setCigar(new Cigar(cigarElements));

        return flipped;
    }
}
