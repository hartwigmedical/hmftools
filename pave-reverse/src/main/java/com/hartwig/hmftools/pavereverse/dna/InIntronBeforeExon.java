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
    public String consistencyWarnings(final GeneTranscript geneTranscript)
    {
        if(!geneTranscript.isCodingBaseAtTheStartOfAnExon(ExonBaseIndex))
        {
            return String.format("Base %d is not at the start of an exon in %s for gene %s",
                    ExonBaseIndex, geneTranscript.Transcript.TransName, geneTranscript.Gene.GeneName);
        }
        return null;
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
