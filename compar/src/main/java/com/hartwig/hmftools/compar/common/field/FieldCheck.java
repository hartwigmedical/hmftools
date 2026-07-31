package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

public class FieldCheck
{
    public final boolean IsCompared;

    public FieldCheck(final boolean isCompared)
    {
        IsCompared = isCompared;
    }

    public String toString() { return format("active(%s)", IsCompared); }
}
