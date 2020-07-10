package com.hartwig.hmftools.svtools.cohort;

public enum LineElementType
{
    KNOWN,
    KNOWN_SUSPECT,
    SUSPECT,
    NONE;

    public static boolean moreKnown(final LineElementType type1, final LineElementType type2)
    {
        if(type1 == type2)
            return false;

        switch(type1)
        {
            case KNOWN:
                return type2 == KNOWN_SUSPECT || type2 == SUSPECT;

            case KNOWN_SUSPECT:
                return type2 == SUSPECT;

            default:
                return false;
        }
    }

    public static LineElementType fromString(final String lineType)
    {
        if(lineType.equals("KNOWN"))
            return LineElementType.KNOWN;
        else if(lineType.contains("KNOWN"))
            return LineElementType.KNOWN_SUSPECT;
        else if(lineType.equals("SUSPECT"))
            return LineElementType.SUSPECT;
        else
            return LineElementType.NONE;
    }

}
