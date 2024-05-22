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

    // indices for the junction sequence around the junction (eg +/- 50 bases)
    private final int mJunctionSeqIndexStart;
    private final int mJunctionSeqIndexEnd;

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
            FullSequence = assembly.formJunctionSequence(RefBaseLength);
        }
        else
        {
            FullSequence = Nucleotides.reverseComplementBases(assembly.formJunctionSequence(RefBaseLength));
        }

        // also make a shorter sequence centred around the junction
        int juncSeqExtLength = maxJuncSeqExtensionLength > 0 ? min(ExtensionLength, maxJuncSeqExtensionLength) : ExtensionLength;
        int juncSeqRefLength = min(RefBaseLength, maxJuncSeqRefBaseLength);

        int juncIndexStart, juncIndexEnd;

        if(assembly.isForwardJunction())
        {
            juncIndexStart = mJunctionIndex - juncSeqRefLength + 1;
            juncIndexEnd = mJunctionIndex + juncSeqExtLength;
        }
        else
        {
            juncIndexStart = mJunctionIndex - juncSeqExtLength;
            juncIndexEnd = mJunctionIndex + juncSeqRefLength - 1;
        }

        if(!Reversed)
        {
            mJunctionSeqIndexStart = juncIndexStart;
            mJunctionSeqIndexEnd = juncIndexEnd;
        }
        else
        {
            // note the switches here
            mJunctionSeqIndexStart = indexReversed(juncIndexEnd);
            mJunctionSeqIndexEnd = indexReversed(juncIndexStart);
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

        if(!Reversed)
        {
            FullSequence = new String(bases);
        }
        else
        {
            FullSequence = Nucleotides.reverseComplementBases(new String(bases));
        }

        // unused
        mJunctionSeqIndexStart = 0;
        mJunctionSeqIndexEnd = 0;
    }

    public int junctionIndex()
    {
        return !Reversed ? mJunctionIndex : BaseLength - mJunctionIndex - 1;
    }

    public final String junctionSequence()
    {
        return FullSequence.substring(junctionSeqStartIndex(), junctionSeqEndIndex() + 1);
    }

    public int junctionSeqStartIndex() { return mJunctionSeqIndexStart; }
    public int junctionSeqEndIndex()
    {
        return mJunctionSeqIndexEnd;
    }

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
        return format("len(%d ref=%d ext=%d juncIndex=%d) %s juncSeq(%d - %d)",
                BaseLength, RefBaseLength, ExtensionLength, mJunctionIndex, Reversed ? "rev" : "fwd",
                junctionSeqStartIndex(), junctionSeqEndIndex());
    }
}
