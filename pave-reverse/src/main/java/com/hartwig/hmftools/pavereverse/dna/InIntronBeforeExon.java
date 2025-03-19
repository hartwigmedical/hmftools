package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public class InIntronBeforeExon implements HgvsAddress
{
    public final int IndexBeforeExonBase;
    public final int ExonBaseIndex;

    public InIntronBeforeExon(int exonBaseIndex, int indexBeforeExonBase)
    {
        IndexBeforeExonBase = indexBeforeExonBase;
        ExonBaseIndex = exonBaseIndex;
    }

    @Override
    public int toStrandLocation(GeneTranscript geneTranscript)
    {
        InExon inExon = new InExon(ExonBaseIndex);
        return inExon.toStrandLocation(geneTranscript) - IndexBeforeExonBase;
    }
}
