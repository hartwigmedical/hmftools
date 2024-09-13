package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import com.hartwig.hmftools.common.genome.region.Orientation;

public class BreakendSegment
{
    public final int AssemblyAlignmentId;
    public final int SequenceIndex;
    public final Orientation Orient;
    public final int Index;

    public final AlignData Alignment;

    public BreakendSegment(
            final int assemblyAlignmentId, final int sequenceIndex, final Orientation orientation, final int index, final AlignData alignment)
    {
        AssemblyAlignmentId = assemblyAlignmentId;
        SequenceIndex = sequenceIndex;
        Orient = orientation;
        Index = index;
        Alignment = alignment;
    }

    public String uniqueId() { return format("%d_%d", AssemblyAlignmentId, Index); }

    public String toString()
    {
        return format("asm(%d index=%d:%d) segment(%d M=%d mq=%d score=%d)",
                AssemblyAlignmentId, SequenceIndex, Orient.asByte(),
                Index, Alignment.alignedBases(), Alignment.mapQual(), Alignment.score());
    }
}
