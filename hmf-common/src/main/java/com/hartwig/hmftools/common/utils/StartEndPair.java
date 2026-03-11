package com.hartwig.hmftools.common.utils;

public enum StartEndPair
{
    Start,
    End;

    public boolean isStart() { return this == Start; }
    public boolean isEnd() { return this == End; }
}
