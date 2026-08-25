package com.hartwig.hmftools.panelbuilder;

// One segment of a probe sequence. Segments are contiguous in probe-sequence space.
// Either a reference genome region or a literal insert sequence.
public sealed interface SequenceSegment extends Comparable<SequenceSegment> permits RefSegment, InsertSeqSegment
{
    // Number of bases this segment contributes to the probe sequence.
    int baseLength();

    // Arbitrary but stable rank so different segment types have a defined order. Distinct per type.
    int typeRank();

    // Compares against another segment of the same concrete type (guaranteed by the equal type rank in compareTo).
    int compareSameType(SequenceSegment other);

    @Override
    default int compareTo(final SequenceSegment other)
    {
        int rank = Integer.compare(typeRank(), other.typeRank());
        return rank != 0 ? rank : compareSameType(other);
    }
}
