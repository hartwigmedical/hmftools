package com.hartwig.hmftools.esvee.sequence;

import java.util.List;

public interface AlignedSequence extends Sequence
{
    List<Alignment> getAlignmentBlocks();
}
