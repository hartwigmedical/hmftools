package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.read.Read;

public class AssemblySupport
{
    private final Read mRead;
    private final SupportType mType;
    private final int mAssemblyIndex;
    private final int mJunctionReadIndex; // index within this read of the junction position
    private final int[] mReadIndexRange;
    private int mJunctionMismatches; // those past the junction
    private int mReferenceMismatches;

    public AssemblySupport(
            final Read read, final SupportType type, final int assemblyIndex, final int junctionReadIndex, final int[] readIndexRange, final int mismatches)
    {
        mRead = read;
        mType = type;
        mAssemblyIndex = assemblyIndex;
        mJunctionReadIndex = junctionReadIndex;
        mReadIndexRange = readIndexRange;
        mJunctionMismatches = mismatches;
        mReferenceMismatches = 0;
    }

    public Read read() { return mRead; }
    public SupportType type() { return mType; }

    public int[] readIndexRange() { return mReadIndexRange; }
    public void setReadIndexRange(final int readIndexStart, final int readIndexEnd)
    {
        mReadIndexRange[SE_START] = readIndexStart;
        mReadIndexRange[SE_END] = readIndexEnd;
    }

    public int readRangeLength() { return mReadIndexRange[1] - mReadIndexRange[0] + 1; }
    public int assemblyIndex() { return mAssemblyIndex; }
    public int junctionReadIndex() { return mJunctionReadIndex; }

    public int mismatchCount() { return mJunctionMismatches + mReferenceMismatches; }
    public int junctionMismatches() { return mJunctionMismatches; }
    public int referenceMismatches() { return mReferenceMismatches; }
    public void setReferenceMismatches(int mismatches) { mReferenceMismatches = mismatches; }

    public static boolean hasMatchingFragment(final List<AssemblySupport> support, final Read read)
    {
        return support.stream().anyMatch(x -> x.read().matchesFragment(read));
    }

    public static List<AssemblySupport> findMatchingFragmentSupport(final List<AssemblySupport> support, final Read read)
    {
        return support.stream().filter(x -> x.read().matchesFragment(read)).collect(Collectors.toList());
    }

    public String toString()
    {
        return format("type(%s) read(%s %d-%d %s) index(asm=%d junc=%d read=%d-%d) mismatch(junc=%d ref=%d)",
                mType, mRead.getName(), mRead.unclippedStart(), mRead.unclippedEnd(), mRead.cigarString(),
                mAssemblyIndex, mJunctionReadIndex, mReadIndexRange[0], mReadIndexRange[1],
                mJunctionMismatches, mReferenceMismatches);
    }
}
