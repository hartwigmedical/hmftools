package com.hartwig.hmftools.esvee.models;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.esvee.processor.SequenceDecomposer;

import htsjdk.samtools.SAMUtils;

public interface Sequence
{
    String getName();

    byte[] getBases();

    default String getBasesString()
    {
        return new String(getBases());
    }

    byte[] getBaseQuality();

    default byte getAverageBaseQuality()
    {
        int sum = 0;
        for(final byte b : getBaseQuality())
            sum += b;
        return (byte) (sum / getLength());
    }

    default int getLength()
    {
        return getBases().length;
    }

    default Sequence subsequence(final int startIndex, final int length)
    {
        return new SimpleSequence(
                Arrays.copyOfRange(getBases(), startIndex, startIndex + length),
                Arrays.copyOfRange(getBaseQuality(), startIndex, startIndex + length));
    }

    default List<SequenceDecomposer.Node> decompose()
    {
        return SequenceDecomposer.decompose(this);
    }

    static Sequence fromBytes(final byte[] bases, final byte[] baseQualities)
    {
        return new SimpleSequence(bases, baseQualities);
    }

    static Sequence fromBytes(final byte[] bases, final byte[] baseQualities, final int length)
    {
        return new SimpleSequence(Arrays.copyOf(bases, length), Arrays.copyOf(baseQualities, length));
    }
}
