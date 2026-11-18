package com.hartwig.hmftools.common.bwa;

import java.util.List;
import java.util.Map;
import java.util.function.Function;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public final class MockBwaMemAligner implements IBwaMemAligner
{
    private final Function<String, List<BwaMemAlignment>> mAlignmentsBySequence;

    // Constructs the aligner such that query sequence maps to a list of alignments.
    // Alignment fails if the query is not in the map.
    public static MockBwaMemAligner forQueries(Map<String, List<BwaMemAlignment>> alignmentsBySequence)
    {
        return new MockBwaMemAligner(sequence ->
        {
            List<BwaMemAlignment> alignments = alignmentsBySequence.get(sequence);
            if(alignments == null)
            {
                throw new IllegalArgumentException("No mock alignments for sequence: " + sequence);
            }
            return alignments;
        });
    }

    // Constructs such that the same list of alignments is returned for all query sequences.
    public static MockBwaMemAligner forAnyQuery(List<BwaMemAlignment> alignments)
    {
        return new MockBwaMemAligner(sequence -> alignments);
    }

    private MockBwaMemAligner(Function<String, List<BwaMemAlignment>> alignmentsBySequence)
    {
        mAlignmentsBySequence = alignmentsBySequence;
    }

    @Override
    public List<List<BwaMemAlignment>> alignSequences(List<byte[]> sequences)
    {
        return sequences.stream().map(sequence -> mAlignmentsBySequence.apply(new String(sequence))).toList();
    }
}
