package com.hartwig.hmftools.pavereverse.util;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Checks
{
    public static final String HGVS_FORMAT_REQUIRED = "Required format is GENE:p.XnY where X and Y are amino acids and n is an integer.";

    public static boolean isNucleotideSequence(String s)
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

    public static boolean isCodon(String s)
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

    public static Matcher matchPattern(final Pattern variationPattern, final String description)
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
