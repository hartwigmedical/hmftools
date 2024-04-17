package com.hartwig.hmftools.sage.common;

import static java.lang.String.format;

public class Microhomology
{
    public final String Bases;
    public final int Length;

    public Microhomology(final String bases, final int length)
    {
        Bases = bases;
        Length = length;
    }

    public String toString()
    {
        return format("%s length(%d)", Bases, Length);
    }
}
