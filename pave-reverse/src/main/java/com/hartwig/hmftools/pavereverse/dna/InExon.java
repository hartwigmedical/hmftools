package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public class InExon implements HgvsAddress
{
    public final int IndexOfBaseInCodingBases;

    public InExon(final int indexOfBaseInCodingBases)
    {
        IndexOfBaseInCodingBases = indexOfBaseInCodingBases;
    }

    @Override
    public int toStrandLocation(GeneTranscript geneTranscript)
    {
        return geneTranscript.absolutePositionOfTranslatedBase(IndexOfBaseInCodingBases);
    }
}
