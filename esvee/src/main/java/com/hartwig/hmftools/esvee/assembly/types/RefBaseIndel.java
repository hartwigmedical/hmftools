package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

public class RefBaseIndel
{
    public final boolean IsInsert;
    public final int Position;
    public final int Length;

    private int mReadSupport;

    public RefBaseIndel(final boolean isInsert, final int position, final int length)
    {
        IsInsert = isInsert;
        Position = position;
        Length = length;
        mReadSupport = 1;
    }

    public int indelLength() { return IsInsert ? Length : -Length; }
    public void addSupport(int support) { mReadSupport += support; }
    public int support() { return mReadSupport; }

    public String toString() { return format("%s %d @ %d", IsInsert ? "insert" : "delete", Length, Position); }
}
