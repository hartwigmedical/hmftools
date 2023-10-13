package com.hartwig.hmftools.common.codon;

public final class Nucleotides
{
    public static final char[] DNA_BASES = {'G', 'A', 'T', 'C'};

    public static char swapDnaBase(final char base)
    {
        if(base == 'A') return 'T';
        if(base == 'T') return 'A';
        if(base == 'C') return 'G';
        if(base == 'G') return 'C';
        return base;
    }

    public static int baseIndex(final char base)
    {
        if(base == 'G')
            return 0;
        if(base == 'A')
            return 1;
        if(base == 'T')
            return 2;
        if(base == 'C')
            return 3;
        else
            return -1;
    }

    public static String swapDnaBase(final String base) { return String.valueOf(swapDnaBase(base.charAt(0))); }

    public static char complement(final char base)
    {
        return swapDnaBase(base);
    }

    public static boolean isValidDnaBase(final char base) { return base == 'G' || base == 'A' || base == 'T' || base == 'C'; }

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
