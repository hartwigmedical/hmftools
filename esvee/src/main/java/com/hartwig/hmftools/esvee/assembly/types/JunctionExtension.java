package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findRepeats;

import java.util.List;

import com.google.common.collect.Lists;

public class JunctionExtension
{
    private final Junction mJunction;

    private byte mBases[];
    private byte mBaseQuals[];
    private final int mRefBaseLength;
    private int mPopulatedBaseLength;

    private final List<SupportRead> mSupport;

    private final List<RepeatInfo> mRepeatInfo;

    private static final int EXTENSION_INCREASE = 200;

    public JunctionExtension(final JunctionAssembly assembly)
    {
        mJunction = assembly.junction();

        mSupport = Lists.newArrayList();
        mRepeatInfo = Lists.newArrayList();

        int extensionLength = assembly.extensionLength();
        mRefBaseLength = min(assembly.refBaseLength(), 30); // need to consider impact of this
        int assemblyCopyLength = extensionLength + mRefBaseLength;
        int initialBaseLength = extensionLength + mRefBaseLength + EXTENSION_INCREASE; // to allow for a few reads to be added

        mBases = new byte[initialBaseLength];
        mBaseQuals = new byte[initialBaseLength];

        // copy ref and extension sequence from the source junction assembly
        if(mJunction.isForward())
        {
            int assemblyCopyIndexStart = assembly.baseLength() - assemblyCopyLength;

            int newBaseIndex = 0;
            for(int assemblyIndex = assemblyCopyIndexStart; assemblyIndex < assembly.baseLength(); ++assemblyIndex, ++newBaseIndex)
            {
                mBases[newBaseIndex] = assembly.bases()[assemblyIndex];
                mBaseQuals[newBaseIndex] = assembly.baseQuals()[assemblyIndex];
            }
        }
        else
        {
            int assemblyCopyIndexEnd = assemblyCopyLength - 1;

            int newBaseIndex = EXTENSION_INCREASE;
            for(int assemblyIndex = 0; assemblyIndex <= assemblyCopyIndexEnd; ++assemblyIndex, ++newBaseIndex)
            {
                mBases[newBaseIndex] = assembly.bases()[assemblyIndex];
                mBaseQuals[newBaseIndex] = assembly.baseQuals()[assemblyIndex];
            }
        }

        mPopulatedBaseLength = assemblyCopyLength;
        buildRepeatInfo();
    }

    public List<RepeatInfo> repeatInfo() { return mRepeatInfo; }

    public void buildRepeatInfo()
    {
        mRepeatInfo.clear();
        List<RepeatInfo> repeats = findRepeats(mBases);
        if(repeats != null)
            mRepeatInfo.addAll(repeats);
    }

    public String toString()
    {
        return format("junc(%s) extensionLen(%d)", mJunction.coords(), mPopulatedBaseLength - mRefBaseLength);
    }

}
