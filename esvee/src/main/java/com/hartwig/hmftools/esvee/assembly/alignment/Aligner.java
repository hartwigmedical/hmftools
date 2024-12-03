package com.hartwig.hmftools.esvee.assembly.alignment;

import java.util.List;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public interface Aligner
{
    List<BwaMemAlignment> alignSequence(final byte[] bases);
    ImmutablePair<List<BwaMemAlignment>,List<BwaMemAlignment>> alignSequences(final byte[] bases1, final byte[] bases2);
}
