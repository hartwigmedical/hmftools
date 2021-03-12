package com.hartwig.hmftools.common.neo;

public enum NeoEpitopeType
{
    MISSENSE,
    INFRAME_INSERTION,
    INFRAME_DELETION,
    FRAMESHIFT,
    STOP_LOST,
    INFRAME_FUSION,
    OUT_OF_FRAME_FUSION;

    public boolean isFusion() { return this == INFRAME_FUSION || this == OUT_OF_FRAME_FUSION; }
    public boolean isPointMutation() { return !isFusion(); }
}
