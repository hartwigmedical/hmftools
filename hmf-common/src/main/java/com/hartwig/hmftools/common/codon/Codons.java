package com.hartwig.hmftools.common.codon;

import org.jetbrains.annotations.NotNull;

public final class Codons
{
    public static final String DNA_START_CODON = "ATG";

    public static final String DNA_STOP_CODON_1 = "TAA";
    public static final String DNA_STOP_CODON_2 = "TAG";
    public static final String DNA_STOP_CODON_3 = "TGA";

    public static final char UNKNOWN = '.';

    public static boolean isStopCodon(final String codon)
    {
        return codon.equals(DNA_STOP_CODON_1) || codon.equals(DNA_STOP_CODON_2) || codon.equals(DNA_STOP_CODON_3);
    }

    public static boolean isStartCodon(final String codon)
    {
        return codon.equals(DNA_START_CODON);
    }

    public static char aminoAcid(final String codon)
    {
        if(isStopCodon(codon))
            return 'X';

        if(isStartCodon(codon))
            return 'M';

        switch(codon)
        {
            // SECOND BASE T
            case "TTT":
            case "TTC":
                return 'F';
            case "TTA":
            case "TTG":
            case "CTT":
            case "CTC":
            case "CTA":
            case "CTG":
                return 'L';
            case "ATT":
            case "ATC":
            case "ATA":
                return 'I';
            case "GTT":
            case "GTC":
            case "GTA":
            case "GTG":
                return 'V';

            // SECOND BASE C
            case "TCT":
            case "TCC":
            case "TCA":
            case "TCG":
                return 'S';
            case "CCT":
            case "CCC":
            case "CCA":
            case "CCG":
                return 'P';
            case "ACT":
            case "ACC":
            case "ACA":
            case "ACG":
                return 'T';
            case "GCT":
            case "GCC":
            case "GCA":
            case "GCG":
                return 'A';

            // SECOND BASE A
            case "TAT":
            case "TAC":
                return 'Y';
            case "CAT":
            case "CAC":
                return 'H';
            case "CAA":
            case "CAG":
                return 'Q';
            case "AAT":
            case "AAC":
                return 'N';
            case "AAA":
            case "AAG":
                return 'K';
            case "GAT":
            case "GAC":
                return 'D';
            case "GAA":
            case "GAG":
                return 'E';

            // SECOND BASE G
            case "TGT":
            case "TGC":
                return 'C';
            case "TGG":
                return 'W';
            case "CGT":
            case "CGC":
            case "CGA":
            case "CGG":
                return 'R';
            case "AGT":
            case "AGC":
                return 'S';
            case "AGA":
            case "AGG":
                return 'R';
            case "GGT":
            case "GGC":
            case "GGA":
            case "GGG":
                return 'G';
        }

        return UNKNOWN;
    }

    public static String codons(@NotNull String aminoAcids)
    {
        StringBuilder builder = new StringBuilder();
        for(int i = 0; i < aminoAcids.length(); i++)
        {
            builder.append(codon(aminoAcids.charAt(i)));
        }
        return builder.toString();
    }

    @NotNull
    public static String codon(char aminoAcid)
    {
        final char[] bases = new char[] { 'G', 'A', 'T', 'C' };
        for(final char firstBase : bases)
        {
            for(final char secondBase : bases)
            {
                for(final char thirdBase : bases)
                {
                    final String codon = String.valueOf(firstBase) + secondBase + thirdBase;
                    if(aminoAcid(codon) == aminoAcid)
                    {
                        return codon;
                    }
                }
            }
        }
        throw new IllegalArgumentException("Unknown amino acid " + aminoAcid);
    }

    @NotNull
    public static String formAminoAcids(@NotNull String dna)
    {
        StringBuilder builder = new StringBuilder();
        for(int i = 0; i < dna.length() - 2; i += 3)
        {
            builder.append(aminoAcid(dna.substring(i, i + 3)));
        }

        return builder.toString();
    }

}
