package com.hartwig.hmftools.common.codon;

public class Nucleotides
{
    public static final char[] DNA_BASES = {'A', 'C', 'G', 'T'};

    public static char swapDnaBase(final char base)
    {
        if(base == 'A') return 'T';
        if(base == 'T') return 'A';
        if(base == 'C') return 'G';
        if(base == 'G') return 'C';
        return base;
    }

    public static String reverseStrandBases(final String bases)
    {
        // reverse and swap base pairs
        String newBases = "";
        for(int i = 0; i < bases.length(); ++i)
        {
            newBases += swapDnaBase(bases.charAt(bases.length() - i - 1));
        }

        return newBases;
    }


}
