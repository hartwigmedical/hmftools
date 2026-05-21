package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.umccr.java.hellbender.utils.bwa.BwaMemAlignment;

public interface PairAligner
{
    ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignSequences(final byte[] bases1, final byte[] bases2);
}
