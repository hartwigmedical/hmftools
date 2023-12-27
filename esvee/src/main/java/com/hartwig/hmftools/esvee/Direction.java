package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;

public enum Direction
{
    FORWARDS(1),
    REVERSE(-1);

    public final int Step;

    Direction(final int step) {
        Step = step;
    }

    public Direction opposite() { return this == FORWARDS ? REVERSE : FORWARDS; }

    public String toShortString() { return this == FORWARDS ? "F" : "R"; }

    public byte toStrand() { return this == FORWARDS ? POS_STRAND : NEG_STRAND; }
}
