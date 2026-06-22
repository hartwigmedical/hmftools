package com.hartwig.hmftools.esvee;

import static java.util.Collections.emptyList;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.bwa.IBwaMemAligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public class MockBwaMemAligner implements IBwaMemAligner
{
    public Map<String, List<BwaMemAlignment>> Alignments;

    public MockBwaMemAligner(final Map<String, List<BwaMemAlignment>> alignments)
    {
        Alignments = alignments;
    }

    @Override
    public List<List<BwaMemAlignment>> alignSequences(final List<byte[]> sequences)
    {
        return sequences.stream().map(this::getAlignments).toList();
    }

    private List<BwaMemAlignment> getAlignments(final byte[] sequence)
    {
        return Alignments.getOrDefault(new String(sequence), emptyList());
    }
}
