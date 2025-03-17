package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public class DownstreamOfCodingEndAddress implements HgvsAddress
{
    public final int IndexDownstreamOfEnd;

    public DownstreamOfCodingEndAddress(final int indexOfBaseInCodingBases)
    {
        this.IndexDownstreamOfEnd = indexOfBaseInCodingBases;
    }

    @Override
    public int toStrandLocation(GeneTranscript geneTranscript)
    {
        return geneTranscript.Transcript.CodingEnd + IndexDownstreamOfEnd;
    }
}
