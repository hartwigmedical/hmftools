package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

public class SequenceDiffInfo
{
    public final int Index;
    public final String Bases; // either the alt SNV base, the indel base or the repeat sequence if a contraction or expansion
    public final SequenceDiffType Type;
    public final int RepeatCount;

    public static final SequenceDiffInfo UNSET = new SequenceDiffInfo(-1, "", SequenceDiffType.UNSET, 0);

    public SequenceDiffInfo(final int index, final String bases, final SequenceDiffType type, final int repeatCount)
    {
        Index = index;
        Bases = bases;
        Type = type;
        RepeatCount = repeatCount;
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
        if(Type == SequenceDiffType.MATCH)
            return format("%d: %s", Index, Type);

        String info = format("%d: %s %s", Index, Type, Bases);
        return RepeatCount == 0 ? info : format("%s x%d)", info, RepeatCount);
    }
}
