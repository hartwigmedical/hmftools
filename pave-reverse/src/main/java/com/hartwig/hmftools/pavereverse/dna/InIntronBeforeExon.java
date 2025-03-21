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
        final int strandLocationOfExonicBase = new InExon(ExonBaseIndex).toStrandLocation(geneTranscript);
        if(geneTranscript.Transcript.posStrand())
        {
            return strandLocationOfExonicBase - IndexBeforeExonBase;
        }
        return strandLocationOfExonicBase + IndexBeforeExonBase;
    }
}
