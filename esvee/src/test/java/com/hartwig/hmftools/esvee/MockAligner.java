package com.hartwig.hmftools.esvee;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.alignment.Aligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public class MockAligner implements Aligner
{
    public final List<List<BwaMemAlignment>> PendingAlignments;

    public MockAligner()
    {
        PendingAlignments = Lists.newArrayList();
    }

    @Override
    public List<BwaMemAlignment> alignSequence(final byte[] bases)
    {
        if(PendingAlignments.isEmpty())
            return Collections.emptyList();

        List<BwaMemAlignment> alignments = PendingAlignments.remove(0);
        return alignments;
    }
}
