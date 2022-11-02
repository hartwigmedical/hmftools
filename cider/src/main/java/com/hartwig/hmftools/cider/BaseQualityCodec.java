package com.hartwig.hmftools.cider;

import static com.google.common.base.Preconditions.checkArgument;

import org.apache.commons.lang3.ArrayUtils;
import org.jetbrains.annotations.NotNull;

// compact way to store base and qualities
// first 2 bit of each byte stores base quality, last 2 bit stores the base
// 6 bit quality can support values up to 64
public class BaseQualityCodec
{
    @NotNull
    public static byte[] encode(@NotNull String readString, @NotNull byte[] baseQualities)
    {
        checkArgument(readString.length() == baseQualities.length);
        byte[] compactDataArray = new byte[readString.length()];
        for (int i = 0; i < readString.length(); ++i)
        {
            compactDataArray[i] = encode(readString.charAt(i), baseQualities[i]);
        }
        return compactDataArray;
    }

    @NotNull
    public static String decodeBaseString(@NotNull byte[] compactDataArray)
    {
        StringBuilder stringBuilder = new StringBuilder(compactDataArray.length);
        for (final byte b : compactDataArray)
        {
            stringBuilder.append(decodeReadBase(b));
        }
        return stringBuilder.toString();
    }

    @NotNull
    public static byte[] decodeBaseQualities(@NotNull byte[] compactDataArray)
    {
        byte[] baseQualities = new byte[compactDataArray.length];
        for (int i = 0; i < compactDataArray.length; ++i)
        {
            baseQualities[i] = decodeBaseQuality(compactDataArray[i]);
        }
        return baseQualities;
    }

    public static char baseAt(@NotNull byte[] compactDataArray, int index)
    {
        return decodeReadBase(compactDataArray[index]);
    }

    public static byte baseQualityAt(@NotNull byte[] compactDataArray, int index)
    {
        return decodeBaseQuality(compactDataArray[index]);
    }

    public static byte encode(char base, byte baseQuality)
    {
        checkArgument(baseQuality <= 63, String.format("illegal baseQuality: %d > 63", baseQuality));

        byte val = (byte) ((baseQuality & 0xFF) << 2);

        switch (base)
        {
            case 'A': break;
            case 'T': val += 3; break;
            case 'C': val += 1; break;
            case 'G': val += 2; break;
            default:
                throw new IllegalArgumentException("illegal base: " + base);
        }

        return val;
    }

    public static char decodeReadBase(byte val)
    {
        switch (val & 0x3)
        {
            case 0: return 'A';
            case 3: return 'T';
            case 1: return 'C';
            case 2: return 'G';
        }

        throw new RuntimeException();
    }

    public static byte decodeBaseQuality(byte val)
    {
        return (byte)((val & 0xFF) >> 2);
    }

    // do the reverse complement inplace
    public static void reverseComplement(byte[] compactDataArray)
    {
        // first reverse the array
        ArrayUtils.reverse(compactDataArray);

        // then we flip each of the base
        for (int i = 0; i < compactDataArray.length; ++i)
        {
            byte data = compactDataArray[i];
            compactDataArray[i] = (byte)((data & 0xFC) + ~data & 0x3);
        }
    }
}
