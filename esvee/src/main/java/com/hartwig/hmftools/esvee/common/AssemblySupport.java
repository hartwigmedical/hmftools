package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import com.hartwig.hmftools.esvee.read.Read;

public class AssemblySupport
{
    private final Read mRead;
    private final int mAssemblyIndex;
    private final int[] mReadIndexRange;
    private final boolean mReferenceBases;
    private final int mMismatches;

    public AssemblySupport(
            final Read read, final int assemblyIndex, final int readIndexStart, final int readIndexEnd, final int mismatches,
            final boolean referenceBases)
    {
        mRead = read;
        mAssemblyIndex = assemblyIndex;
        mReadIndexRange = new int[] { readIndexStart, readIndexEnd };
        mMismatches = mismatches;
        mReferenceBases = referenceBases;
    }

    public Read read() { return mRead; }
    public int mismatches() { return mMismatches; }
    public int[] readIndexRange() { return mReadIndexRange; }
    public int readRangeLength() { return mReadIndexRange[1] - mReadIndexRange[0] + 1; }
    public int assemblyIndex() { return mAssemblyIndex; }
    public boolean isReferenceBases() { return mReferenceBases; }

    public String toString()
    {
        return format("read(%s) assIndex(%d) readIndeRange(%d-%d) mismatch(%d) type(%s)",
                mRead.getName(), mAssemblyIndex, mReadIndexRange[0], mReadIndexRange[1], mMismatches,
                mReferenceBases ? "ref" : "non-ref");
    }
}
