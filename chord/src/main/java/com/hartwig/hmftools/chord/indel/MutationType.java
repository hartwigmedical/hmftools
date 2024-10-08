package com.hartwig.hmftools.chord.indel;

public enum MutationType
{
    DELETION("del"),
    INSERTION("ins");

    public final String Abbreviation;

    MutationType(String abbreviation)
    {
        Abbreviation = abbreviation;
    }
}
