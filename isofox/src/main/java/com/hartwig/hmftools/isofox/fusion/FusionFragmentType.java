package com.hartwig.hmftools.isofox.fusion;

public enum FusionFragmentType
{
    UNKNOWN,
    MATCHED_JUNCTION,
    REALIGN_CANDIDATE,
    REALIGNED,
    DISCORDANT,
    DISCORDANT_JUNCTION;

    public boolean isJunctionType() { return this == MATCHED_JUNCTION || this == DISCORDANT_JUNCTION; }
}
