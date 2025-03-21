package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public class InIntronAfterExon implements HgvsAddress
{
    public final int IndexAfterExonBase;
    public final int ExonBaseIndex;

    public InIntronAfterExon(int exonBaseIndex, int indexAfterExonBase)
    {
        IndexAfterExonBase = indexAfterExonBase;
        ExonBaseIndex = exonBaseIndex;
    }

    @Override
    public int toStrandLocation(GeneTranscript geneTranscript)
    {
        final int strandLocationOfExonicBase = new InExon(ExonBaseIndex).toStrandLocation(geneTranscript);
        if(geneTranscript.Transcript.posStrand())
        {
            return strandLocationOfExonicBase + IndexAfterExonBase;
        }
        return strandLocationOfExonicBase - IndexAfterExonBase;
    }
}
