package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import java.util.Set;

public class SbxDuplexIndel
{
    public final int DuplexIndelIndexStart;
    public final int DuplexIndelIndexEnd;
    public final String RepeatBases;
    public final int FirstReadInsertIndex;

    public final Set<Integer> LowBaseQualIndices; // after base trimming

    public final int DeletedIndelIndexStart;
    public final int DeletedIndelIndexEnd;

    public SbxDuplexIndel(
            int duplexIndelIndexStart, int duplexIndelIndexEnd, String repeatBases,
            int firstReadInsertIndex, final Set<Integer> lowBaseQualIndices, int deletedIndelIndexStart, int deletedIndelIndexEnd)
    {
        DuplexIndelIndexStart = duplexIndelIndexStart;
        DuplexIndelIndexEnd = duplexIndelIndexEnd;
        RepeatBases = repeatBases;
        FirstReadInsertIndex = firstReadInsertIndex;
        LowBaseQualIndices = lowBaseQualIndices;
        DeletedIndelIndexStart = deletedIndelIndexStart;
        DeletedIndelIndexEnd = deletedIndelIndexEnd;
    }

    public int duplexLowQualCount() { return DuplexIndelIndexEnd - DuplexIndelIndexStart + 1; }
    public int deletedBaseCount() { return DeletedIndelIndexEnd - DeletedIndelIndexStart + 1; }

    public boolean withinBounds(int readIndex)
    {
        return readIndex >= FirstReadInsertIndex && readIndex <= DuplexIndelIndexEnd;
    }

    public boolean withinDeleteBounds(int readIndex)
    {
        return readIndex >= DeletedIndelIndexStart && readIndex <= DeletedIndelIndexEnd;
    }

    public boolean isLowQualBase(int readIndex)
    {
        return LowBaseQualIndices.contains(readIndex);
    }

    public String toString()
    {
        return format("duplexIndel(%d-%d len=%d) repeat(%s start=%d) lowBq(%d) deleted(%d-%d)",
                DuplexIndelIndexStart, DuplexIndelIndexEnd, duplexLowQualCount(), RepeatBases,
                FirstReadInsertIndex, LowBaseQualIndices.size(), DeletedIndelIndexStart, DeletedIndelIndexEnd);
    }
}
