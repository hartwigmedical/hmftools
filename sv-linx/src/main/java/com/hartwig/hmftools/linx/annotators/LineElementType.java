package com.hartwig.hmftools.linx.annotators;

import java.util.Set;

public enum LineElementType
{
    NONE,
    KNOWN,
    SUSPECT;

    public static String toString(final Set<LineElementType> types)
    {
        if(types.isEmpty())
            return NONE.toString();
        else if(types.size() == 1)
            return types.iterator().next().toString();
        else
            return String.format("%s;%s", KNOWN.toString(), SUSPECT.toString());

    }
}
