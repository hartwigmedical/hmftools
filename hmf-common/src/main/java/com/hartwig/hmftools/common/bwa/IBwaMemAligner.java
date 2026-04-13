package com.hartwig.hmftools.common.bwa;

import java.util.List;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public interface IBwaMemAligner
{
    List<List<BwaMemAlignment>> alignSequences(final List<byte[]> sequences);

    default List<BwaMemAlignment> alignSequence(final byte[] sequence)
    {
        return alignSequences(List.of(sequence)).get(0);
    }
}
