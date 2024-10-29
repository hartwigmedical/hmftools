package com.hartwig.hmftools.esvee.assembly.alignment;

import static java.lang.String.format;

public class BreakendSegment
{
    public final int AssemblyAlignmentId;
    public final int SequenceIndex; // the index in the full alignment sequence where this segment starts
    public final int Index; // the index of this sequence out of all segments in the full assembly

    public final AlignData Alignment;

    private final int[] mIndelSeqenceIndices; // if from a CIGAR indel, then the full sequence indices of the indel

    public BreakendSegment(
            final int assemblyAlignmentId, final int sequenceIndex, final int index, final AlignData alignment)
    {
        this(assemblyAlignmentId, sequenceIndex, index, alignment, null);
    }

    public BreakendSegment(
            final int assemblyAlignmentId, final int sequenceIndex, final int index, final AlignData alignment, final int[] indelIndices)
    {
        AssemblyAlignmentId = assemblyAlignmentId;
        SequenceIndex = sequenceIndex;
        Index = index;
        Alignment = alignment;
        mIndelSeqenceIndices = indelIndices;
    }

    public int[] indelSeqenceIndices() { return mIndelSeqenceIndices; }

    public String uniqueId() { return format("%d_%d", AssemblyAlignmentId, Index); }

    public String toString()
    {
        return format("%d: seqIndex(%d) alignment(M=%d mq=%d score=%d)",
                Index, SequenceIndex, Alignment.alignedBases(), Alignment.mapQual(), Alignment.score());
    }
}
