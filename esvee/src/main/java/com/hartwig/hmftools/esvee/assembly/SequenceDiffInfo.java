package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.String.format;

public class SequenceDiffInfo
{
    public final int Index;
    public final SequenceDiffType Type;
    public final int BaseCount; // positive means was an extra base in the first sequence and vice versa

    public SequenceDiffInfo(final int index, final SequenceDiffType type, final int baseCount)
    {
        Index = index;
        Type = type;
        BaseCount = baseCount;
    }

    public boolean isFirst() { return BaseCount > 0; }

    public String toString() { return format("%d: %s %s len(%d)", Index, Type, isFirst() ? "1st" : "2nd", abs(BaseCount)); }
}
