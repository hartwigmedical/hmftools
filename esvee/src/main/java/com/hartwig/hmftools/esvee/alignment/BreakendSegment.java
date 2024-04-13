package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

public class BreakendSegment
{
    public final int AssemblyAlignmentId;
    public final int SequenceLength;
    public final int SequenceIndex;
    public final byte Orientation;
    public final int Index;

    public final AlignData Alignment;

    public BreakendSegment(
            final int assemblyAlignmentId, final int sequenceLength, final int sequenceIndex, final byte orientation,
            final int index, final AlignData alignment)
    {
        AssemblyAlignmentId = assemblyAlignmentId;
        SequenceLength = sequenceLength;
        SequenceIndex = sequenceIndex;
        Orientation = orientation;
        Index = index;
        Alignment = alignment;
    }

    public String toString()
    {
        return format("asm(%d len=%d index=%d:%d) segment(%d M=%d mq=%d score=%d trimLen=%d)",
                AssemblyAlignmentId, SequenceLength, SequenceIndex, Orientation,
                Index, Alignment.alignedBases(), Alignment.MapQual, Alignment.Score, Alignment.repeatTrimmedLength());
    }
}
