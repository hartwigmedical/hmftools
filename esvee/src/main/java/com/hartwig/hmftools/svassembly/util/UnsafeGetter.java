package com.hartwig.hmftools.svassembly.util;

import java.lang.reflect.Field;

import sun.misc.Unsafe;

public enum UnsafeGetter
{
    ;

    public static final Unsafe UNSAFE;

    static
    {
        final Field f;
        try
        {
            f = Unsafe.class.getDeclaredField("theUnsafe");
            f.setAccessible(true);
            UNSAFE = (Unsafe) f.get(null);
        }
        catch(final Exception e)
        {
            throw new RuntimeException(e);
        }
    }
}
