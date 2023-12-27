package com.hartwig.hmftools.esvee.util;

public final class CommonUtils
{
    // could be moved into BaseRegion as common methods
    public static boolean touches(final int leftStart, final int leftEnd, final int rightStart, final int rightEnd)
    {
        return leftStart - 1 <= rightEnd && leftEnd + 1 >= rightStart;
    }

    public static boolean overlaps(final int leftStart, final int leftEnd, final int rightStart, final int rightEnd)
    {
        return leftStart <= rightEnd && leftEnd >= rightStart;
    }

    public static byte[] reverseBytes(final byte[] bytes)
    {
        if(bytes == null)
            return null;

        final byte[] reversed = new byte[bytes.length];
        for(int i = 0; i < bytes.length; i++)
        {
            reversed[reversed.length - i - 1] = bytes[i];
        }

        return reversed;
    }

    public static String formatNanos(final long nanoseconds)
    {
        final String[] unitNames = new String[] { "ns", "us", "ms", "sec", "min", "hr", "days" };
        final int[] divisors = new int[] { 1000, 1000, 1000, 60, 60, 24 };

        int unitIndex = 0;
        double time = nanoseconds;

        while(unitIndex < divisors.length && time > 5 * divisors[unitIndex])
            time /= divisors[unitIndex++];

        return String.format("%.2f %s", time, unitNames[unitIndex]);
    }
}
