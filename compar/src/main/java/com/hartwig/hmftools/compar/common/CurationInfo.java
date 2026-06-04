package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

public class CurationInfo
{
    public final CurationType Type;
    public final String Comment;

    public CurationInfo(final CurationType type, final String comment)
    {
        Type = type;
        Comment = comment;
    }

    public String toString() { return format("%s: %s", Type, Comment); }
}
