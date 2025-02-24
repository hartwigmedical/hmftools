package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.common.codon.AminoAcids.AMINO_ACID_TO_CODON_MAP;

import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.AminoAcids;

import org.jetbrains.annotations.NotNull;

class AminoAcid
{
    static final AminoAcid START = new AminoAcid("M");

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

    @NotNull
    public Set<String> matchingCodons(@NotNull final String prefix, @NotNull final String suffix)
    {
        Preconditions.checkArgument(prefix.length() < 3);
        Preconditions.checkArgument(suffix.length() < 3);
        return matchingCodons(prefix, suffix, false);
    }

    @NotNull
    public Set<String> matchingTruncatedCodons(@NotNull final String prefix, @NotNull final String suffix)
    {
        Preconditions.checkArgument(prefix.length() < 3);
        Preconditions.checkArgument(suffix.length() < 3);
        return matchingCodons(prefix, suffix, true);
    }

    @NotNull
    private Set<String> matchingCodons(@NotNull final String prefix, @NotNull final String suffix, final boolean truncate)
    {
        final Set<String> result = new HashSet<>();
        AMINO_ACID_TO_CODON_MAP.get(symbol).forEach(codon ->
        {
            if(codon.startsWith(prefix) && codon.endsWith(suffix))
            {
                if(truncate)
                {
                    String prefixTruncatedCodon = codon.substring(prefix.length());
                    String truncatedCodon = prefixTruncatedCodon.substring(0, prefixTruncatedCodon.length() - suffix.length());
                    if(!truncatedCodon.isEmpty())
                    {
                        result.add(truncatedCodon);
                    }
                }
                else
                {
                    result.add(codon);
                }
            }
        });
        return result;
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
