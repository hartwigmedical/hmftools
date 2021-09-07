package com.hartwig.hmftools.common.codon;

import java.util.List;

import com.google.common.collect.Lists;

public class Nucleotides
{
    public static final char[] DNA_BASES = {'G', 'A', 'T', 'C'};
    public static final List<String> DNA_BASE_LIST = Lists.newArrayList("G", "A", "T", "C");

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
        StringBuilder newBases = new StringBuilder();

        for(int i = bases.length() - 1; i >= 0; i--)
        {
            newBases.append(swapDnaBase(bases.charAt(i)));
        }

        return newBases.toString();
    }

}
