package com.hartwig.hmftools.common.genome.region;

import org.jetbrains.annotations.NotNull;

public enum Strand
{
    FORWARD,
    REVERSE;

    // alternative byte representation for file I/O
    public static final byte POS_STRAND = 1;
    public static final byte NEG_STRAND = -1;

    @NotNull
    public static Strand valueOf(int direction)
    {
        switch(direction)
        {
            case 1:
                return Strand.FORWARD;
            case -1:
                return Strand.REVERSE;
        }

        throw new IllegalArgumentException("Invalid direction: " + direction);
    }

    public static Strand fromChar(char strand)
    {
        switch(strand)
        {
            case '+':
                return Strand.FORWARD;
            case '-':
                return Strand.REVERSE;
        }

        throw new IllegalArgumentException("Invalid strand: " + strand);
    }

    public char asChar() { return this == FORWARD ? '+' : '-'; }

    public Strand getOpposite()
    {
        switch(this)
        {
            case FORWARD:
                return Strand.REVERSE;
            case REVERSE:
                return Strand.FORWARD;
        }
        throw new IllegalArgumentException("Invalid strand: " + this);
    }
}
