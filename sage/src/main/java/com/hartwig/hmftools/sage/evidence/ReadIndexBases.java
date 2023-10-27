package com.hartwig.hmftools.sage.evidence;

public class ReadIndexBases
{
    public final int Index;
    public final byte[] Bases;

    public ReadIndexBases(final int index, final byte[] bases)
    {
        Index = index;
        Bases = bases;
    }
}
