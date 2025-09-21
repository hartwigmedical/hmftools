package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

public class SbxDuplexIndel
{
    public final int DuplexIndelIndexStart;
    public final int DuplexIndelIndexEnd;
    public final String RepeatBases;
    public final int ReadRepeatIndexStart;
    public final int RefRepeatCount;
    public final int DeletedIndelIndexStart;
    public final int DeletedIndelIndexEnd;

    public SbxDuplexIndel(
            final int duplexIndelIndexStart, final int duplexIndelIndexEnd, final String repeatBases, final int refRepeatCount,
            final int readRepeatIndexStart, final int deletedIndelIndexStart, final int deletedIndelIndexEnd)
    {
        DuplexIndelIndexStart = duplexIndelIndexStart;
        DuplexIndelIndexEnd = duplexIndelIndexEnd;
        RepeatBases = repeatBases;
        ReadRepeatIndexStart = readRepeatIndexStart;
        RefRepeatCount = refRepeatCount;
        DeletedIndelIndexStart = deletedIndelIndexStart;
        DeletedIndelIndexEnd = deletedIndelIndexEnd;
    }

    public int deletedBaseCount() { return DeletedIndelIndexEnd - DeletedIndelIndexStart + 1; }

    public String toString()
    {
        return format("duplexIndel(%d-%d) repeat(%s index=%d) refCount(%d) deleted(%d-%d)",
                DuplexIndelIndexStart, DuplexIndelIndexEnd, RepeatBases, ReadRepeatIndexStart, RefRepeatCount,
                DeletedIndelIndexStart, DeletedIndelIndexEnd);
    }
}
