package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.AminoAcids;

import org.jetbrains.annotations.NotNull;

class AminoAcid
{
    @NotNull
    final String symbol;

    static boolean isValidAminoAcidIdentifier(String s)
    {
        if(AminoAcids.AMINO_ACID_TO_CODON_MAP.containsKey(s))
        {
            return true;
        }
        return AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.containsKey(s);
    }


    public AminoAcid(@NotNull final String symbol)
    {
        Preconditions.checkArgument(isValidAminoAcidIdentifier(symbol), "Not a valid amino acid identifier: " + symbol);
        this.symbol = AminoAcids.forceSingleLetterProteinAnnotation(symbol);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final AminoAcid aminoAcid = (AminoAcid) o;
        return Objects.equals(symbol, aminoAcid.symbol);
    }

    @Override
    public int hashCode()
    {
        return Objects.hashCode(symbol);
    }

    @Override
    public String toString()
    {
        return symbol;
    }
}
