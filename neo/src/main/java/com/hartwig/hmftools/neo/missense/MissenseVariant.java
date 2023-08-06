package com.hartwig.hmftools.neo.missense;

import static java.lang.String.format;

public class MissenseVariant
{
    public String TransName;
    public int Position; // of the missense variant
    public int CodonIndex;
    public String Context;
    public char RefBase;
    public char AltBase;

    public MissenseVariant(
            final String transName, final int position, final int codonIndex, final String context, final char refBase, final char altBase)
    {
        TransName = transName;
        Position = position;
        CodonIndex = codonIndex;
        Context = context;
        RefBase = refBase;
        AltBase = altBase;
    }

    public String toString()
    {
        return format("%s var(%d %c>%c) codon(%d)", Position, RefBase, AltBase, CodonIndex);
    }

}
