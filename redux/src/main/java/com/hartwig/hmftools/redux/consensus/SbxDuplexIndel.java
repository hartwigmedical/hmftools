package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

public class SbxDuplexIndel
{
    public final int DuplexIndelIndexStart;
    public final int DuplexIndelIndexEnd;
    public final String RepeatBases;
    public final int TotalRepeatBaseLength;
    public final int FirstReadInsertIndex;

    public final List<Integer> LowBaseQualIndices; // after base trimming

    public final int DeletedIndelIndexStart;
    public final int DeletedIndelIndexEnd;

    public SbxDuplexIndel(
            int duplexIndelIndexStart, int duplexIndelIndexEnd, String repeatBases, int totalRepeatBaseLength,
            int firstReadInsertIndex, final List<Integer> lowBaseQualIndices, int deletedIndelIndexStart, int deletedIndelIndexEnd)
    {
        DuplexIndelIndexStart = duplexIndelIndexStart;
        DuplexIndelIndexEnd = duplexIndelIndexEnd;
        RepeatBases = repeatBases;
        TotalRepeatBaseLength = totalRepeatBaseLength;
        FirstReadInsertIndex = firstReadInsertIndex;
        LowBaseQualIndices = lowBaseQualIndices;
        DeletedIndelIndexStart = deletedIndelIndexStart;
        DeletedIndelIndexEnd = deletedIndelIndexEnd;
    }

    public SbxDuplexIndel(final SbxDuplexIndel other, int readOffset)
    {
        DuplexIndelIndexStart = other.DuplexIndelIndexStart + readOffset;
        DuplexIndelIndexEnd = other.DuplexIndelIndexEnd + readOffset;
        RepeatBases = other.RepeatBases;
        TotalRepeatBaseLength = other.TotalRepeatBaseLength;
        FirstReadInsertIndex = other.FirstReadInsertIndex + readOffset;
        LowBaseQualIndices = Lists.newArrayListWithCapacity(other.LowBaseQualIndices.size());
        other.LowBaseQualIndices.forEach(x -> LowBaseQualIndices.add(x + readOffset));
        DeletedIndelIndexStart = other.DeletedIndelIndexStart + readOffset;
        DeletedIndelIndexEnd = other.DeletedIndelIndexEnd + readOffset;
    }

    public int duplexLowQualCount() { return DuplexIndelIndexEnd - DuplexIndelIndexStart + 1; }
    public int deletedBaseCount() { return DeletedIndelIndexEnd - DeletedIndelIndexStart + 1; }

    public boolean isLowQualBase(int readIndex)
    {
        return LowBaseQualIndices.contains(readIndex);
    }

    public String toString()
    {
        return format("duplexIndel(%d-%d len=%d) repeat(%s full-len=%d start=%d) lowBq(%d) deleted(%d-%d)",
                DuplexIndelIndexStart, DuplexIndelIndexEnd, duplexLowQualCount(), RepeatBases, TotalRepeatBaseLength,
                FirstReadInsertIndex, LowBaseQualIndices.size(), DeletedIndelIndexStart, DeletedIndelIndexEnd);
    }
}
