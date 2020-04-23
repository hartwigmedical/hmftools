package com.hartwig.hmftools.isofox.fusion;

public enum FusionFragmentType
{
    UNKNOWN,
    BOTH_JUNCTIONS,
    REALIGNED,
    ONE_JUNCTION, // may resolve to one of the others
    ONE_SIDED,
    DISCORDANT;
}
