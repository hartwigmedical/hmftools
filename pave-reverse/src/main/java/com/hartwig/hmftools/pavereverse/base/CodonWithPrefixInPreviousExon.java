package com.hartwig.hmftools.pavereverse.base;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;

class CodonWithPrefixInPreviousExon extends CodonWithinExons
{
    private final String mPrefix;

    CodonWithPrefixInPreviousExon(String prefix, BaseSequence body)
    {
        super(body);
        mPrefix = prefix;
    }

    public int forwardStrandLocationOfChange(int positionWithinCodonOfChange)
    {
        Preconditions.checkArgument(positionWithinCodonOfChange >= mPrefix.length());
        Preconditions.checkArgument(positionWithinCodonOfChange <= 2);
        return Body.Start + positionWithinCodonOfChange - mPrefix.length();
    }

    @Override
    public CodonWithinExons reverseComplement()
    {
        return new CodonWithSuffixInNextExon(Body.reverseComplement(), Nucleotides.reverseComplementBases(mPrefix));
    }

    @Override
    public String fixedPrefix()
    {
        return mPrefix;
    }
}
