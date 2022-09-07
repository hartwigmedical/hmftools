package com.hartwig.hmftools.common.codon;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASES;

import org.jetbrains.annotations.NotNull;

public final class Codons
{
    public static final String START_CODON = "ATG";

    public static final String STOP_CODON_1 = "TAA";
    public static final String STOP_CODON_2 = "TAG";
    public static final String STOP_CODON_3 = "TGA";

    public static final char START_AMINO_ACID = 'M';
    public static final char STOP_AMINO_ACID = 'X';

    public static final char UNKNOWN = '.';

    private static final char[] sCodonLookup = {'K', 'N', 'N', 'K', 'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'R', 'S', 'S', 'R',
            'X', 'Y', 'Y', 'X', 'L', 'F', 'F', 'L', 'S', 'S', 'S', 'S', 'X', 'C', 'C', 'W', 'Q', 'H', 'H', 'Q', 'L', 'L', 'L', 'L',
            'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'E', 'D', 'D', 'E', 'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G'};

    public static boolean isCodonMultiple(int bases) { return (bases % 3) == 0; }

    public static boolean isStopCodon(final String codon)
    {
        return codon.equals(STOP_CODON_1) || codon.equals(STOP_CODON_2) || codon.equals(STOP_CODON_3);
    }

    public static boolean isStartCodon(final String codon)
    {
        return codon.equals(START_CODON);
    }

    public static char codonToAminoAcid(final String codon)
    {
        return codonToAminoAcid(codon, 0);
    }

    public static char codonToAminoAcid(final String dna, int index)
    {
        try
        {
            // convert the 3 bases into a look up index
            int lookupIndex = baseToInt(dna.charAt(index)) * 16 +
                    baseToInt(dna.charAt(index + 1)) * 4 +
                    baseToInt(dna.charAt(index + 2));
            if (lookupIndex < sCodonLookup.length)
                return sCodonLookup[lookupIndex];
        }
        catch (IllegalArgumentException ignored) {}
        return UNKNOWN;
    }

    public static String aminoAcidsToCodons(@NotNull String aminoAcids)
    {
        StringBuilder builder = new StringBuilder();
        for(int i = 0; i < aminoAcids.length(); i++)
        {
            builder.append(aminoAcidToCodon(aminoAcids.charAt(i)));
        }
        return builder.toString();
    }

    @NotNull
    public static String aminoAcidToCodon(char aminoAcid)
    {
        for(final char firstBase : DNA_BASES)
        {
            for(final char secondBase : DNA_BASES)
            {
                for(final char thirdBase : DNA_BASES)
                {
                    final String codon = String.valueOf(firstBase) + secondBase + thirdBase;
                    if(codonToAminoAcid(codon) == aminoAcid)
                    {
                        return codon;
                    }
                }
            }
        }

        throw new IllegalArgumentException("Unknown amino acid " + aminoAcid);
    }

    // note: this is used only for the lookup table. The ordering must be
    // kept this way otherwise the lookup table needs to be recreated
    private static int baseToInt(char base)
    {
        switch (base)
        {
            case 'A':
                return 0;
            case 'T':
                return 1;
            case 'C':
                return 2;
            case 'G':
                return 3;
        }
        throw new IllegalArgumentException("unknown base: " + base);
    }

    @NotNull
    public static String aminoAcidFromBases(@NotNull String dna)
    {
        StringBuilder builder = new StringBuilder();
        for(int i = 0; i < dna.length() - 2; i += 3)
        {
            builder.append(codonToAminoAcid(dna, i));
        }

        return builder.toString();
    }
}
