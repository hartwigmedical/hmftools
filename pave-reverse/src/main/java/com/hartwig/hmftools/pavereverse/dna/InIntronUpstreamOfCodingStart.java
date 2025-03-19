package com.hartwig.hmftools.pavereverse.dna;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public class InIntronUpstreamOfCodingStart implements HgvsAddress
{
    public final int IndexOfExonicBase;
    public final int RelativePositionOfIntronicBase;

    public InIntronUpstreamOfCodingStart(int indexOfExonicBase, final int relativePositionOfIntronicBase)
    {
        Preconditions.checkArgument(indexOfExonicBase < 0);
        RelativePositionOfIntronicBase = relativePositionOfIntronicBase;
        IndexOfExonicBase = indexOfExonicBase;
    }

    @Override
    public int toStrandLocation(GeneTranscript geneTranscript)
    {
        return new InExonUpstreamOfCodingStart(IndexOfExonicBase).toStrandLocation(geneTranscript) + RelativePositionOfIntronicBase;
    }
}
