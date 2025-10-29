package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public class InExonDownstreamOfCodingEnd implements HgvsAddress
{
    public final int IndexDownstreamOfEnd;

    public InExonDownstreamOfCodingEnd(final int indexOfBaseInCodingBases)
    {
        IndexDownstreamOfEnd = indexOfBaseInCodingBases;
    }

    @Override
    public int toStrandLocation(GeneTranscript geneTranscript)
    {
        return geneTranscript.absolutePositionOf3PrimeUtrExonicBase(IndexDownstreamOfEnd);
    }
}
