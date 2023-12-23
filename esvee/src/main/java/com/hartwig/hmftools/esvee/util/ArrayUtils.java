package com.hartwig.hmftools.esvee.util;

import org.jetbrains.annotations.Nullable;

public enum ArrayUtils
{
    ;

    @Nullable
    public static byte[] reverse(final byte[] bytes) {
        if (bytes == null)
            return null;

        final byte[] reversed = new byte[bytes.length];
        for (int i = 0; i < bytes.length; i++)
            reversed[reversed.length - i - 1] = bytes[i];
        return reversed;
    }
}
