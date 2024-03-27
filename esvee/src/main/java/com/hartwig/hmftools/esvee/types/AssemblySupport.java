package com.hartwig.hmftools.esvee.types;

import static java.lang.String.format;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.read.Read;

public class AssemblySupport
{
    private final Read mRead;
    private final SupportType mType;
    private final int mAssemblyIndex;
    private final int mJunctionReadIndex; // index within this read of the junction position
    // private final int[] mReadIndexRange;

    // those past the junction
    private int mJunctionMatches;
    private int mJunctionMismatches;
    private int mReferenceMismatches;

    public AssemblySupport(
            final Read read, final SupportType type, final int assemblyIndex, final int junctionReadIndex,
            final int matches, final int mismatches)
    {
        mRead = read;
        mType = type;
        mAssemblyIndex = assemblyIndex;
        mJunctionReadIndex = junctionReadIndex;
        mJunctionMatches = matches;
        mJunctionMismatches = mismatches;
        mReferenceMismatches = 0;
    }

    public Read read() { return mRead; }
    public SupportType type() { return mType; }

    public int assemblyIndex() { return mAssemblyIndex; }
    public int junctionReadIndex() { return mJunctionReadIndex; }

    public int mismatchCount() { return mJunctionMismatches + mReferenceMismatches; }
    public int junctionMismatches() { return mJunctionMismatches; }
    public int junctionMatches() { return mJunctionMatches; }
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
        return format("type(%s) read(%s %d-%d %s) index(asm=%d junc=%d) hqMatch(%d) mismatch(junc=%d ref=%d)",
                mType, mRead.getName(), mRead.unclippedStart(), mRead.unclippedEnd(), mRead.cigarString(),
                mAssemblyIndex, mJunctionReadIndex, mJunctionMatches, mJunctionMismatches, mReferenceMismatches);
    }

    /*
    // public int readRangeLength() { return mReadIndexRange[1] - mReadIndexRange[0] + 1; }
    public int[] readIndexRange() { return mReadIndexRange; }
    public void setReadIndexRange(final int readIndexStart, final int readIndexEnd)
    {
        mReadIndexRange[SE_START] = readIndexStart;
        mReadIndexRange[SE_END] = readIndexEnd;
    }
    */
}
