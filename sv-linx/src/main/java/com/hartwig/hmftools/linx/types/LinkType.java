package com.hartwig.hmftools.linx.types;

public enum LinkType
{
    TEMPLATED_INSERTION,
    DELETION_BRIDGE;

    public static String linkTypeStr(final LinkType type)
    {
        return type == TEMPLATED_INSERTION ? "TI" : "DB";
    }
}
