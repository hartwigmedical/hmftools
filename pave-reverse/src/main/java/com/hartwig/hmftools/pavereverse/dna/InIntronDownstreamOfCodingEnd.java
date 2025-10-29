package com.hartwig.hmftools.pavereverse.dna;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public class InIntronDownstreamOfCodingEnd implements HgvsAddress
{
    public final int IndexOfExonicBase;
    public final int RelativePositionOfIntronicBase;

    public InIntronDownstreamOfCodingEnd(int indexOfExonicBase, int relativePositionOfIntronicBase)
    {
        Preconditions.checkArgument(indexOfExonicBase > 0);
        RelativePositionOfIntronicBase = relativePositionOfIntronicBase;
        IndexOfExonicBase = indexOfExonicBase;
    }

    @Override
    public int toStrandLocation(GeneTranscript geneTranscript)
    {
        final int exonicBaseLocation = new InExonDownstreamOfCodingEnd(IndexOfExonicBase).toStrandLocation(geneTranscript);
        if(geneTranscript.Transcript.posStrand())
        {
            return exonicBaseLocation + RelativePositionOfIntronicBase;
        }
        return exonicBaseLocation - RelativePositionOfIntronicBase;
    }
}
