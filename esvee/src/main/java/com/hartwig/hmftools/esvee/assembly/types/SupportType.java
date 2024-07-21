package com.hartwig.hmftools.esvee.assembly.types;

public enum SupportType
{
    JUNCTION,
    INDEL,
    DISCORDANT,
    JUNCTION_MATE,
    EXTENSION;

    public boolean isSplitSupport() { return this == JUNCTION || this == INDEL; }
}
