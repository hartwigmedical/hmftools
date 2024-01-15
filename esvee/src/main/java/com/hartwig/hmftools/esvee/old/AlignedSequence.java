package com.hartwig.hmftools.esvee.old;

import java.util.List;

public interface AlignedSequence extends Sequence
{
    List<Alignment> getAlignmentBlocks();
}
