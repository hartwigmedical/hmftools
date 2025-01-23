package com.hartwig.hmftools.pave.transval;

import com.hartwig.hmftools.common.codon.AminoAcids;

import org.jetbrains.annotations.NotNull;

public class Checks
{
    static boolean isNucleotideSequence(@NotNull String s)
    {
        if(s.isEmpty())
        {
            return false;
        }
        for(int i = 0; i < s.length(); i++)
        {
            if(!isNucleotide(s.charAt(0)))
            {
                return false;
            }
        }
        return true;
    }

    static boolean isValidAminoAcidLetter(String s)
    {
        return AminoAcids.AMINO_ACID_TO_CODON_MAP.containsKey(s);
    }

    static boolean isCodon(@NotNull String s)
    {
        if(s.length() != 3)
        {
            return false;
        }
        return isNucleotide(s.charAt(0)) && isNucleotide(s.charAt(1)) && isNucleotide(s.charAt(2));
    }

    static boolean isNucleotide(char c)
    {
        return c == 'A' || c == 'C' || c == 'G' || c == 'T';
    }

    static boolean isValidAminoAcidName(String s)
    {
        if(AminoAcids.AMINO_ACID_TO_CODON_MAP.containsKey(s))
        {
            return true;
        }
        return (AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.containsKey(s));
    }
}
