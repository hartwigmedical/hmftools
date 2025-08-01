package com.hartwig.hmftools.common.codon;

public final class Nucleotides
{
    public static final char[] DNA_BASES = {'A', 'C', 'G', 'T'};
    public static final byte[] DNA_BASE_BYTES = { 65, 67, 71, 84 };
    public static final char DNA_N_BASE = 'N';
    public static final byte DNA_N_BYTE = 78;

    public static char swapDnaBase(final char base)
    {
        if(base == 'A') return 'T';
        if(base == 'T') return 'A';
        if(base == 'C') return 'G';
        if(base == 'G') return 'C';
        return base;
    }

    public static byte swapDnaBase(final byte base) { return (byte)swapDnaBase((char)base); }

    public static int baseIndex(final byte base)
    {
        if(base == DNA_BASE_BYTES[0]) return 0;
        if(base == DNA_BASE_BYTES[1]) return 1;
        if(base == DNA_BASE_BYTES[2]) return 2;
        if(base == DNA_BASE_BYTES[3]) return 3;
        if(base == DNA_N_BYTE) return 4;

        return -1;
    }

    public static int baseIndex(final char base)
    {
        if(base == DNA_BASES[0]) return 0;
        if(base == DNA_BASES[1]) return 1;
        if(base == DNA_BASES[2]) return 2;
        if(base == DNA_BASES[3]) return 3;
        if(base == DNA_N_BASE) return 4;

        return -1;
    }

    public static String swapDnaBase(final String base) { return String.valueOf(swapDnaBase(base.charAt(0))); }

    public static char complement(final char base)
    {
        return swapDnaBase(base);
    }

    public static boolean isValidDnaBase(final char base)
    {
        return base == DNA_BASES[0] || base == DNA_BASES[1] || base == DNA_BASES[2] || base == DNA_BASES[3];
    }

    public static boolean isValidDnaBase(final byte base)
    {
        return base == DNA_BASE_BYTES[0] || base == DNA_BASE_BYTES[1] || base == DNA_BASE_BYTES[2] || base == DNA_BASE_BYTES[3];
    }

    public static String reverseComplementBases(final String bases)
    {
        // reverse and swap base pairs
        StringBuilder newBases = new StringBuilder();

        for(int i = bases.length() - 1; i >= 0; i--)
        {
            newBases.append(swapDnaBase(bases.charAt(i)));
        }

        return newBases.toString();
    }

    public static byte[] reverseComplementBases(final byte[] bases)
    {
        byte[] newBases = new byte[bases.length];

        for(int i = 0; i < bases.length; ++i)
        {
            newBases[i] = swapDnaBase(bases[bases.length - 1 - i]);
        }

        return newBases;
    }
}
