package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public class UpstreamOfCodingStartAddress implements HgvsAddress
{
    public final int IndexUpstreamOfStart;

    public UpstreamOfCodingStartAddress(final int indexOfBaseInCodingBases)
    {
        this.IndexUpstreamOfStart = indexOfBaseInCodingBases;
    }

    @Override
    public int toStrandLocation(GeneTranscript geneTranscript)
    {
        return geneTranscript.Transcript.CodingStart - IndexUpstreamOfStart;
    }
}
