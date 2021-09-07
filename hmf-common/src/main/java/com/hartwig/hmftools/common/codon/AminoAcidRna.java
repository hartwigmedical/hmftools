package com.hartwig.hmftools.common.codon;

import static com.hartwig.hmftools.common.codon.Codons.STOP_CODON_1;
import static com.hartwig.hmftools.common.codon.Codons.STOP_CODON_2;
import static com.hartwig.hmftools.common.codon.Codons.STOP_CODON_3;
import static com.hartwig.hmftools.common.codon.Codons.codonToAminoAcid;
import static com.hartwig.hmftools.common.codon.Codons.isStopCodon;

public class AminoAcidRna
{
    private static final String RNA_STOP_CODON_1 = swapDnaToRna(STOP_CODON_1);
    private static final String RNA_STOP_CODON_2 = swapDnaToRna(STOP_CODON_2);
    private static final String RNA_STOP_CODON_3 = swapDnaToRna(STOP_CODON_3);

    public static final String STOP_SYMBOL = "_";

    public static final String START_CODON = "AUG";
    public static final String AA_SELENOCYSTEINE = "U";


    public static boolean isRnaStopCodon(final String rnaCodon)
    {
        return rnaCodon.equals(RNA_STOP_CODON_1) || rnaCodon.equals(RNA_STOP_CODON_2) || rnaCodon.equals(RNA_STOP_CODON_3);
    }

    public static char swapRnaBase(final char base)
    {
        if(base == 'A') return 'U';
        if(base == 'U') return 'A';
        if(base == 'C') return 'G';
        if(base == 'G') return 'C';
        return base;
    }

    public static String swapDnaToRna(final String bases)
    {
        return bases.replaceAll("T", "U");
    }

    public static String swapRnaToDna(final String bases)
    {
        return bases.replaceAll("U", "T");
    }

    public static String convertDnaCodonToAminoAcid(final String codon)
    {
        if(isStopCodon(codon))
            return STOP_SYMBOL;

        return String.valueOf(codonToAminoAcid(codon));
    }

    public static String convertRnaCodonToAminoAcid(final String codon)
    {
        if(isRnaStopCodon(codon))
            return STOP_SYMBOL;

        return String.valueOf(codonToAminoAcid(swapRnaToDna(codon)));
    }
}
