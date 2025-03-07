package com.hartwig.hmftools.pavereverse.base;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;

import org.jetbrains.annotations.NotNull;

class CodonWithPrefixInPreviousExon extends CodonWithinExons
{
    @NotNull
    private final String Prefix;

    CodonWithPrefixInPreviousExon(@NotNull String prefix, final BaseSequence body)
    {
        super(body);
        this.Prefix = prefix;
    }

    public int forwardStrandLocationOfChange(int positionWithinCodonOfChange)
    {
        Preconditions.checkArgument(positionWithinCodonOfChange >= Prefix.length());
        Preconditions.checkArgument(positionWithinCodonOfChange <= 2);
        return Body.Start + positionWithinCodonOfChange - Prefix.length();
    }

    @Override
    @NotNull
    public CodonWithinExons reverseComplement()
    {
        return new CodonWithSuffixInNextExon(Body.reverseComplement(), Nucleotides.reverseComplementBases(Prefix));
    }

    @NotNull
    @Override
    public String fixedPrefix()
    {
        return Prefix;
    }
}
