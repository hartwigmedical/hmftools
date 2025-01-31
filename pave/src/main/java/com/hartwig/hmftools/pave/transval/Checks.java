package com.hartwig.hmftools.pave.transval;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.hartwig.hmftools.common.codon.AminoAcids;

import org.jetbrains.annotations.NotNull;

public class Checks
{
    static final String HGVS_FORMAT_REQUIRED = "Required format is GENE:p.XnY where X and Y are amino acids and n is an integer.";

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

    static boolean isValidProtein(String s)
    {
        if(s.isEmpty())
        {
            return false;
        }
        for(int i = 0; i < s.length(); i++)
        {
            if (!isValidAminoAcidLetter("" + s.charAt(i)))
            {
                return false;
            }
        }
        return true;
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

    @NotNull
    static Matcher matchPattern(final Pattern variationPattern, final String description)
    {
        final Matcher matcher = variationPattern.matcher(description);
        boolean matches = matcher.find();
        if(!matches)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        return matcher;
    }
}
