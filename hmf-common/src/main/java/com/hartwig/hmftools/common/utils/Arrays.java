package com.hartwig.hmftools.common.utils;

import static java.lang.System.arraycopy;

import static com.hartwig.hmftools.common.sigs.DataUtils.doublesEqual;

import com.hartwig.hmftools.common.codon.Nucleotides;

public final class Arrays
{
    public static byte[] subsetArray(final byte[] source, final int startIndex, final int endIndex)
    {
        int length = endIndex - startIndex + 1;
        byte[] dest = new byte[length];
        arraycopy(source, startIndex, dest, 0, length);
        return dest;
    }

    public static void copyArray(final byte[] source, final byte[] dest, final int sourceIndexStart, final int destIndexStart)
    {
        int d = destIndexStart;
        for(int s = sourceIndexStart; s < source.length && d < dest.length; ++s, ++d)
        {
            dest[d] = source[s];
        }
    }

    public static void copyArray(
            final byte[] source, final byte[] dest, final int sourceIndexStart, final int sourceIndexEnd, final int destIndexStart)
    {
        int d = destIndexStart;
        for(int s = sourceIndexStart; s < source.length & s < sourceIndexEnd && d < dest.length; ++s, ++d)
        {
            dest[d] = source[s];
        }
    }

    public static byte[] copyArray(final byte[] source)
    {
        byte[] dest = new byte[source.length];
        arraycopy(source, 0, dest, 0, source.length);
        return dest;
    }

    public static int[] copyArray(final int[] source)
    {
        int[] dest = new int[source.length];
        arraycopy(source, 0, dest, 0, source.length);
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

    public static boolean equalArray(double[] array1, double[] array2)
    {
        if(array1.length != array2.length)
            return false;

        for(int i = 0; i < array1.length; ++i)
        {
            if(!doublesEqual(array1[i], array2[i]))
                return false;
        }

        return true;
    }

    public static boolean equalArray(byte[] array1, byte[] array2)
    {
        if(array1.length != array2.length)
            return false;

        for(int i = 0; i < array1.length; ++i)
        {
            if(array1[i] != array2[i])
                return false;
        }

        return true;
    }

    public static boolean equalArray(int[] array1, int[] array2)
    {
        if(array1.length != array2.length)
            return false;

        for(int i = 0; i < array1.length; ++i)
        {
            if(array1[i] != array2[i])
                return false;
        }

        return true;
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
