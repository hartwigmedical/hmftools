package com.hartwig.hmftools.pavereverse.base;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;

import org.jetbrains.annotations.NotNull;

class CodonWithSuffixInNextExon extends CodonWithinExons
{
    @NotNull
    private final String Suffix;

    CodonWithSuffixInNextExon(final BaseSequence body, @NotNull String suffix)
    {
        super(body);
        this.Suffix = suffix;
    }

    public int forwardStrandLocationOfChange(int positionWithinCodonOfChange)
    {
        Preconditions.checkArgument(positionWithinCodonOfChange >= 0);
        Preconditions.checkArgument(positionWithinCodonOfChange <= 2 - Suffix.length());
        return Body.Start + positionWithinCodonOfChange;
    }

    @NotNull
    @Override
    public String fixedSuffix()
    {
        return Suffix;
    }

    @Override
    @NotNull
    public CodonWithinExons reverseComplement()
    {
        return new CodonWithPrefixInPreviousExon(Nucleotides.reverseComplementBases(Suffix), Body.reverseComplement());
    }
}
