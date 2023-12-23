package com.hartwig.hmftools.svassembly.util;

public enum StringUtils
{
    ;

    public static String toSnakeCase(final String input)
    {
        final StringBuilder sb = new StringBuilder(input.length() * 2);
        boolean lastWasLower = false;
        int upperRun = 0;
        for(int i = 0; i < input.length(); i++)
        {
            final char c = input.charAt(i);
            if(Character.isUpperCase(c))
            {
                if(lastWasLower)
                {
                    sb.append("_");
                    lastWasLower = false;
                }
                upperRun++;
            }
            else
            {
                if(upperRun > 1)
                {
                    final char lastChar = sb.charAt(sb.length() - 1);
                    sb.deleteCharAt(sb.length() - 1);
                    sb.append("_");
                    sb.append(lastChar);
                }
                upperRun = 0;
                lastWasLower = true;
            }
            sb.append(Character.toLowerCase(c));
        }

        return sb.toString();
    }

    public static boolean parseBoolean(final String value)
    {
        final String lowerValue = value.toLowerCase();
        if(lowerValue.equals("true") || lowerValue.equals("t") || lowerValue.equals("1"))
            return true;
        else if(lowerValue.equals("false") || lowerValue.equals("f") || lowerValue.equals("0") || lowerValue.equals("-1"))
            return false;

        throw new IllegalArgumentException("Cannot parse value as boolean: " + value);
    }

    public static String toSubscriptString(final int n)
    {
        final String string = String.valueOf(n);
        final StringBuilder sb = new StringBuilder(string.length());
        for(int i = 0; i < string.length(); i++)
        {
            final int digit = string.charAt(i) - '0';
            sb.append((char) ('₀' + digit));
        }
        return sb.toString();
    }

    public static String formatNanos(final long nanoseconds)
    {
        final String[] unitNames = new String[] { "ns", "us", "ms", "sec", "min", "hr", "days" };
        final int[] divisors = new int[] { 1000, 1000, 1000, 60, 60, 24 };

        int unitIndex = 0;
        double time = nanoseconds;

        while (unitIndex < divisors.length && time > 5 * divisors[unitIndex])
            time /= divisors[unitIndex++];

        return String.format("%.2f %s", time, unitNames[unitIndex]);
    }
}
