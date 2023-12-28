package com.hartwig.hmftools.esvee.models;

import java.util.List;

public interface AlignedSequence extends Sequence
{
    List<Alignment> getAlignmentBlocks();
}
