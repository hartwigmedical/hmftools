package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import com.hartwig.hmftools.esvee.read.Read;

public class AssemblySupport
{
    private final Read mRead;
    private final int mAssemblyIndex;
    private final int[] mReadIndexRange;
    private final int mMismatches;

    public AssemblySupport(final Read read, final int assemblyIndex, final int readIndexStart, final int readIndexEnd, final int mismatches)
    {
        mRead = read;
        mAssemblyIndex = assemblyIndex;
        mReadIndexRange = new int[] { readIndexStart, readIndexEnd };
        mMismatches = mismatches;
    }

    public String toString()
    {
        return format("read(%s) assIndex(%d) readIndeRange(%d-%d) mismatch(%d)",
                mRead.getName(), mAssemblyIndex, mReadIndexRange[0], mReadIndexRange[1], mMismatches);
    }
}
