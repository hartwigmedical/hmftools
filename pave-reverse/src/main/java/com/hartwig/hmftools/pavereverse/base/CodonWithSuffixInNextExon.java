package com.hartwig.hmftools.pavereverse.base;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;

class CodonWithSuffixInNextExon extends CodonWithinExons
{
    private final String Suffix;

    CodonWithSuffixInNextExon(final BaseSequence body, String suffix)
    {
        super(body);
        Suffix = suffix;
    }

    public int forwardStrandLocationOfChange(int positionWithinCodonOfChange)
    {
        Preconditions.checkArgument(positionWithinCodonOfChange >= 0);
        Preconditions.checkArgument(positionWithinCodonOfChange <= 2 - Suffix.length());
        return Body.Start + positionWithinCodonOfChange;
    }

    @Override
    public String fixedSuffix()
    {
        return Suffix;
    }

    @Override
    public CodonWithinExons reverseComplement()
    {
        return new CodonWithPrefixInPreviousExon(Nucleotides.reverseComplementBases(Suffix), Body.reverseComplement());
    }
}
