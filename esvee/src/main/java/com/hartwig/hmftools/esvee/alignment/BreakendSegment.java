package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import com.hartwig.hmftools.common.genome.region.Orientation;

public class BreakendSegment
{
    public final int AssemblyAlignmentId;
    public final int SequenceLength;
    public final int SequenceIndex;
    public final Orientation Orient;
    public final int Index;

    public final AlignData Alignment;

    public BreakendSegment(
            final int assemblyAlignmentId, final int sequenceLength, final int sequenceIndex, final Orientation orientation,
            final int index, final AlignData alignment)
    {
        AssemblyAlignmentId = assemblyAlignmentId;
        SequenceLength = sequenceLength;
        SequenceIndex = sequenceIndex;
        Orient = orientation;
        Index = index;
        Alignment = alignment;
    }

    public int calcQual()
    {
        int repeatAdjustment = Alignment.segmentLength() - Alignment.repeatTrimmedLength();
        double lengthFactor = (Alignment.Score - repeatAdjustment)/100.0;
        return (int)round(Alignment.MapQual * min(1, lengthFactor));
    }

    public String toString()
    {
        return format("asm(%d len=%d index=%d:%d) segment(%d M=%d mq=%d score=%d trimLen=%d)",
                AssemblyAlignmentId, SequenceLength, SequenceIndex, Orient,
                Index, Alignment.alignedBases(), Alignment.MapQual, Alignment.Score, Alignment.repeatTrimmedLength());
    }
}
