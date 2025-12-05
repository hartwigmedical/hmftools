package com.hartwig.hmftools.esvee.assembly;

public enum BaseQualType
{
    HIGH,
    MEDIUM,
    LOW;

    public static BaseQualType selectHigher(final BaseQualType first, final BaseQualType second)
    {
        return first.ordinal() <= second.ordinal() ? first : second;
    }

    public static BaseQualType selectLower(final BaseQualType first, final BaseQualType second)
    {
        return first.ordinal() >= second.ordinal() ? first : second;
    }
}
