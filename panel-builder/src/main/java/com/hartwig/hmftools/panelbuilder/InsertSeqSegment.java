package com.hartwig.hmftools.panelbuilder;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.SequenceUtils.isDnaSequenceNormal;

// A probe segment consisting of a literal inserted nucleotide sequence (i.e. not based on the reference genome).
public record InsertSeqSegment(String sequence) implements SequenceSegment
{
    public InsertSeqSegment
    {
        requireNonNull(sequence);
        if(sequence.isEmpty())
        {
            throw new IllegalArgumentException("Empty insert sequence");
        }
        if(!isDnaSequenceNormal(sequence))
        {
            throw new IllegalArgumentException("Invalid insert sequence: " + sequence);
        }
    }

    @Override
    public int baseLength()
    {
        return sequence.length();
    }

    @Override
    public int typeRank()
    {
        return 1;
    }

    @Override
    public int compareSameType(final SequenceSegment other)
    {
        return sequence.compareTo(((InsertSeqSegment) other).sequence);
    }
}
