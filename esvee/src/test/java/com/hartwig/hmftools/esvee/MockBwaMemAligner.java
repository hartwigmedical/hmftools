package com.hartwig.hmftools.esvee;

import static java.util.Collections.emptyList;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.bwa.IBwaMemAligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public class MockBwaMemAligner implements IBwaMemAligner
{
    // Use empty string for alignment output for all query sequences.
    public Map<String, List<BwaMemAlignment>> Alignments;

    public MockBwaMemAligner(final Map<String, List<BwaMemAlignment>> alignments)
    {
        Alignments = alignments;
    }

    public MockBwaMemAligner()
    {
        this(new HashMap<>());
    }

    @Override
    public List<List<BwaMemAlignment>> alignSequences(final List<byte[]> sequences)
    {
        return sequences.stream().map(this::getAlignments).toList();
    }

    private List<BwaMemAlignment> getAlignments(final byte[] sequence)
    {
        List<BwaMemAlignment> sequenceAlignments = Alignments.getOrDefault(new String(sequence), emptyList());
        List<BwaMemAlignment> genericAlignments = Alignments.getOrDefault("", emptyList());
        return Stream.concat(sequenceAlignments.stream(), genericAlignments.stream()).toList();
    }
}
