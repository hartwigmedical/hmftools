package com.hartwig.hmftools.esvee.models;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;

public interface AlignedSequence extends Sequence
{
    List<Alignment> getAlignmentBlocks();
}
