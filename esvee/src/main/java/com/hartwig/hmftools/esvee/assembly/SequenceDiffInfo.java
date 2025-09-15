package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

public class SequenceDiffInfo
{
    public final int Index;
    public final String Bases; // either the alt SNV base, the indel base or the repeat sequence if a contraction or expansion
    public final SequenceDiffType Type;
    public final int RepeatDiff;

    public SequenceDiffInfo(final int index, final String bases, final SequenceDiffType type, final int repeatDiff)
    {
        Index = index;
        Bases = bases;
        Type = type;
        RepeatDiff = repeatDiff;
    }

    public static SequenceDiffInfo fromSnv(final int index, final byte base)
    {
        return new SequenceDiffInfo(index, String.valueOf((char)base), SequenceDiffType.BASE, 0);
    }

    public static SequenceDiffInfo fromMatch(final int index)
    {
        return new SequenceDiffInfo(index, null, SequenceDiffType.MATCH, (short)0);
    }

    public String toString()
    {
        String info = format("%d: %s %s", Index, Bases, Type);
        return RepeatDiff == 0 ? info : format("%s repeatDiff(%d)", info, RepeatDiff);
    }
}
