package com.hartwig.hmftools.esvee;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.esvee.alignment.Aligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public class MockAligner implements Aligner
{
    public MockAligner()
    {
    }

    @Override
    public List<BwaMemAlignment> alignSequence(final byte[] bases)
    {
        return Collections.emptyList();
    }

}
