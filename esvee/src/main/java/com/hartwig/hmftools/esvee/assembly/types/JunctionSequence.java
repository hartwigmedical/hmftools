package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.AssemblyConstants.PHASED_ASSEMBLY_JUNCTION_OVERLAP;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class JunctionSequence
{
    public final boolean Reversed;
    public final String FullSequence;

    public final int ExtensionLength;
    public final int RefBaseLength; // may be capped
    public final int BaseLength;

    private byte[] mOriginalBases;
    private byte[] mOriginalBaseQuals;

    // the following are in non-reversed terms
    private final int mJunctionIndex; // as per the original assembly

    // indices for the sequence used to match to another, typically around the junction (eg +/- 50 bases)
    private final int mMatchSeqIndexStart;
    private final int mMatchSeqIndexEnd;

    // built on demand since only used for the sequence comparison routine
    private List<RepeatInfo> mRepeatInfo;
    private byte[] mBases;
    private byte[] mBaseQuals;

    public JunctionSequence(final JunctionAssembly assembly, final boolean reverseCompliment)
    {
        this(assembly, reverseCompliment, PHASED_ASSEMBLY_JUNCTION_OVERLAP, PHASED_ASSEMBLY_JUNCTION_OVERLAP);
    }

    public JunctionSequence(
            final JunctionAssembly assembly, final boolean reverseCompliment, final int maxJuncSeqRefBaseLength, final int maxJuncSeqExtensionLength)
    {
        mOriginalBases = assembly.bases();
        mOriginalBaseQuals = assembly.baseQuals();

        // these have lazy initialisation and aren't used if there is a simple (string search) sequence match between assemblies
        mBases = null;
        mBaseQuals = null;
        mRepeatInfo = null;

        Reversed = reverseCompliment;
        RefBaseLength = assembly.refBaseLength();
        ExtensionLength = assembly.extensionLength();
        BaseLength = RefBaseLength + ExtensionLength;

        mJunctionIndex = assembly.junctionIndex();

        if(!Reversed)
        {
            FullSequence = assembly.formFullSequence();
        }
        else
        {
            FullSequence = Nucleotides.reverseComplementBases(assembly.formFullSequence());
        }

        // also make a shorter sequence centred around the junction
        int matchSeqExtLength = maxJuncSeqExtensionLength > 0 ? min(ExtensionLength, maxJuncSeqExtensionLength) : ExtensionLength;
        int matchSeqRefLength = min(RefBaseLength, maxJuncSeqRefBaseLength);

        int matchIndexStart, matchIndexEnd;

        if(assembly.isForwardJunction())
        {
            matchIndexStart = mJunctionIndex - matchSeqRefLength + 1;
            matchIndexEnd = mJunctionIndex + matchSeqExtLength;
        }
        else
        {
            matchIndexStart = mJunctionIndex - matchSeqExtLength;
            matchIndexEnd = mJunctionIndex + matchSeqRefLength - 1;
        }

        if(!Reversed)
        {
            mMatchSeqIndexStart = matchIndexStart;
            mMatchSeqIndexEnd = matchIndexEnd;
        }
        else
        {
            // note the switches here
            mMatchSeqIndexStart = indexReversed(matchIndexEnd);
            mMatchSeqIndexEnd = indexReversed(matchIndexStart);
        }
    }

    public JunctionSequence(final byte[] bases, final byte[] baseQuals, final Orientation orientation,final boolean reverseCompliment)
    {
        mOriginalBases = bases;
        mOriginalBaseQuals = baseQuals;

        // these have lazy initialisation and aren't used if there is a simple (string search) sequence match between assemblies
        mBases = null;
        mBaseQuals = null;
        mRepeatInfo = null;

        Reversed = reverseCompliment;
        RefBaseLength = 0;
        ExtensionLength = bases.length;
        BaseLength = RefBaseLength + ExtensionLength;

        mJunctionIndex = orientation.isReverse() ? 0 : bases.length - 1; // better to not set this??

        FullSequence = !Reversed ? new String(bases) : new String(Nucleotides.reverseComplementBases(bases));

        // unused
        mMatchSeqIndexStart = 0;
        mMatchSeqIndexEnd = 0;
    }

    public int junctionIndex()
    {
        return !Reversed ? mJunctionIndex : BaseLength - mJunctionIndex - 1;
    }

    public final String matchSequence()
    {
        return FullSequence.substring(mMatchSeqIndexStart, mMatchSeqIndexEnd + 1);
    }

    public int matchSeqStartIndex() { return mMatchSeqIndexStart; }
    public int matchSeqEndIndex() { return mMatchSeqIndexEnd; }

    public byte[] originalBases() { return mOriginalBases; }
    public byte[] originalBaseQuals() { return mOriginalBaseQuals; }

    public byte[] bases()
    {
        if(!Reversed)
        {
            return mOriginalBases;
        }

        if(mBases == null)
        {
            mBases = FullSequence.getBytes();
        }

        return mBases;
    }

    public byte[] baseQuals()
    {
        if(!Reversed)
        {
            return mOriginalBaseQuals;
        }

        if(mBaseQuals == null)
        {
            int baseLength = mOriginalBaseQuals.length;
            mBaseQuals = new byte[baseLength];

            for(int i = 0; i < baseLength; ++i)
            {
                mBaseQuals[i] = mOriginalBaseQuals[baseLength - i - 1];
            }
        }

        return mBaseQuals;
    }

    public List<RepeatInfo> repeatInfo()
    {
        if(mRepeatInfo == null)
        {
            List<RepeatInfo> repeats = RepeatInfo.findRepeats(FullSequence.getBytes());
            mRepeatInfo = repeats != null ? repeats : Collections.emptyList();
        }

        return mRepeatInfo;
    }

    private int indexReversed(int index)
    {
        int juncIndexDiff = index - mJunctionIndex;
        return junctionIndex() - juncIndexDiff;
    }

    public int indexReverted(int index)
    {
        int juncIndexDiff = index - junctionIndex();
        return mJunctionIndex - juncIndexDiff;
    }

    public String toString()
    {
        return format("len(%d ref=%d ext=%d juncIndex=%d) %s matchSeq(%d - %d)",
                BaseLength, RefBaseLength, ExtensionLength, mJunctionIndex, Reversed ? "rev" : "fwd",
                matchSeqStartIndex(), matchSeqEndIndex());
    }
}
