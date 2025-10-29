package com.hartwig.hmftools.pavereverse.dna;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public class InExonUpstreamOfCodingStart implements HgvsAddress
{
    public final int IndexUpstreamOfStart;

    public InExonUpstreamOfCodingStart(final int indexOfBaseInCodingBases)
    {
        Preconditions.checkArgument(indexOfBaseInCodingBases < 0);
        IndexUpstreamOfStart = indexOfBaseInCodingBases;
    }

    @Override
    public int toStrandLocation(GeneTranscript geneTranscript)
    {
        return geneTranscript.absolutePositionOf5PrimeUtrExonicBase(IndexUpstreamOfStart);
    }
}
