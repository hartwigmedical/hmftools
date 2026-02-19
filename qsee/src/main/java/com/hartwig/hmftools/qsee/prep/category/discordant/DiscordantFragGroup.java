package com.hartwig.hmftools.qsee.prep.category.discordant;

import com.hartwig.hmftools.common.sv.DiscordantFragType;

public enum DiscordantFragGroup
{
    DEL_DUP_SHORT("Del/Dup 1K-5K"),
    DEL_DUP_MEDIUM("Del/Dup 5K-100K"),
    DEL_DUP_LONG("Del/Dup >100K"),
    INV_SHORT("Inv <=1K"),
    INV_MEDIUM("Inv 1K-100K"),
    INV_LONG("Inv >100K"),
    TRANSLOCATION("Translocation");

    public final String mName;

    DiscordantFragGroup(String name)
    {
        mName = name;
    }

    public String getName(){ return mName; }

    public static DiscordantFragGroup fromType(DiscordantFragType type)
    {
        return switch(type)
        {
            case Del1To5K, Dup1To5K -> DiscordantFragGroup.DEL_DUP_SHORT;
            case Del5To100K, Dup5To100K -> DiscordantFragGroup.DEL_DUP_MEDIUM;
            case DelGt100K, DupGt100K -> DiscordantFragGroup.DEL_DUP_LONG;
            case InvLt1K -> DiscordantFragGroup.INV_SHORT;
            case Inv1To5K, Inv5To100K -> DiscordantFragGroup.INV_MEDIUM;
            case InvGt100K -> DiscordantFragGroup.INV_LONG;
            case Translocation -> DiscordantFragGroup.TRANSLOCATION;
        };
    }
}
