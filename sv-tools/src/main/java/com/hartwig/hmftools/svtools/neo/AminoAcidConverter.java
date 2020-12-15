package com.hartwig.hmftools.svtools.neo;

public class AminoAcidConverter
{
    public static final String STOP_CODON_1 = "UAA";
    public static final String STOP_CODON_2 = "UAG";
    public static final String STOP_CODON_3 = "UGA";

    public static final String STOP_SYMBOL = "_";

    public static final String START_CODON = "AUG";

    public static final String UNKNOWN = "ERROR";

    public static boolean isStopCodon(final String codon)
    {
        final String rnaCodon = swapDnaToRna(codon);

        return rnaCodon.equals(STOP_CODON_1) || rnaCodon.equals(STOP_CODON_2) || rnaCodon.equals(STOP_CODON_3);
    }

    public static char swapDnaBase(final char base)
    {
        if(base == 'A') return 'T';
        if(base == 'T') return 'A';
        if(base == 'C') return 'G';
        if(base == 'G') return 'C';
        return base;
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

    public static String convertDnaCodonToAminoAcid(final String dnaCodon)
    {
        return convertCodonToAminoAcid(swapDnaToRna(dnaCodon));
    }

    public static String convertCodonToAminoAcid(final String codon)
    {
        if(isStopCodon(codon))
            return STOP_SYMBOL;

        if(codon.equals(START_CODON)) return "M";

        if(codon.equals("UCA")) return "S";
        if(codon.equals("AUA")) return "I";
        if(codon.equals("UCC")) return "S";
        if(codon.equals("AUC")) return "I";
        if(codon.equals("UCG")) return "S";
        if(codon.equals("AUU")) return "I";
        if(codon.equals("UCU")) return "S";
        if(codon.equals("UUC")) return "F";
        if(codon.equals("ACA")) return "T";
        if(codon.equals("UUU")) return "F";
        if(codon.equals("ACC")) return "T";
        if(codon.equals("UUA")) return "L";
        if(codon.equals("ACG")) return "T";
        if(codon.equals("UUG")) return "L";
        if(codon.equals("ACU")) return "T";
        if(codon.equals("UAC")) return "Y";
        if(codon.equals("AAC")) return "N";
        if(codon.equals("UAU")) return "Y";
        if(codon.equals("AAU")) return "N";
        if(codon.equals("AAA")) return "K";
        if(codon.equals("AAG")) return "K";
        if(codon.equals("UGC")) return "C";
        if(codon.equals("AGC")) return "S";
        if(codon.equals("UGU")) return "C";
        if(codon.equals("AGU")) return "S";
        if(codon.equals("AGA")) return "R";
        if(codon.equals("UGG")) return "W";
        if(codon.equals("AGG")) return "R";
        if(codon.equals("CUA")) return "L";
        if(codon.equals("GUA")) return "V";
        if(codon.equals("CUC")) return "L";
        if(codon.equals("GUC")) return "V";
        if(codon.equals("CUG")) return "L";
        if(codon.equals("GUG")) return "V";
        if(codon.equals("CUU")) return "L";
        if(codon.equals("GUU")) return "V";
        if(codon.equals("CCA")) return "P";
        if(codon.equals("GCA")) return "A";
        if(codon.equals("CCC")) return "P";
        if(codon.equals("GCC")) return "A";
        if(codon.equals("CCG")) return "P";
        if(codon.equals("GCG")) return "A";
        if(codon.equals("CCU")) return "P";
        if(codon.equals("GCU")) return "A";
        if(codon.equals("CAC")) return "H";
        if(codon.equals("GAC")) return "D";
        if(codon.equals("CAU")) return "H";
        if(codon.equals("GAU")) return "D";
        if(codon.equals("CAA")) return "Q";
        if(codon.equals("GAA")) return "E";
        if(codon.equals("CAG")) return "Q";
        if(codon.equals("GAG")) return "E";
        if(codon.equals("CGA")) return "R";
        if(codon.equals("GGA")) return "G";
        if(codon.equals("CGC")) return "R";
        if(codon.equals("GGC")) return "G";
        if(codon.equals("CGG")) return "R";
        if(codon.equals("GGG")) return "G";
        if(codon.equals("CGU")) return "R";
        if(codon.equals("GGU")) return "G";

        return UNKNOWN;
    }

    public static String convertAminoAcidToDnaCodon(final String aminoAcid)
    {
        return swapRnaToDna(convertAminoAcidToCodon(aminoAcid));
    }

    public static String convertAminoAcidToCodon(final String aminoAcid)
    {
        if(aminoAcid.equals("S")) return "UCA";
        if(aminoAcid.equals("I")) return "AUA";
        if(aminoAcid.equals("S")) return "UCC";
        if(aminoAcid.equals("I")) return "AUC";
        if(aminoAcid.equals("S")) return "UCG";
        if(aminoAcid.equals("I")) return "AUU";
        if(aminoAcid.equals("S")) return "UCU";
        if(aminoAcid.equals("M")) return "AUG";
        if(aminoAcid.equals("F")) return "UUC";
        if(aminoAcid.equals("T")) return "ACA";
        if(aminoAcid.equals("F")) return "UUU";
        if(aminoAcid.equals("T")) return "ACC";
        if(aminoAcid.equals("L")) return "UUA";
        if(aminoAcid.equals("T")) return "ACG";
        if(aminoAcid.equals("L")) return "UUG";
        if(aminoAcid.equals("T")) return "ACU";
        if(aminoAcid.equals("Y")) return "UAC";
        if(aminoAcid.equals("N")) return "AAC";
        if(aminoAcid.equals("Y")) return "UAU";
        if(aminoAcid.equals("N")) return "AAU";
        if(aminoAcid.equals("K")) return "AAA";
        if(aminoAcid.equals("K")) return "AAG";
        if(aminoAcid.equals("C")) return "UGC";
        if(aminoAcid.equals("S")) return "AGC";
        if(aminoAcid.equals("C")) return "UGU";
        if(aminoAcid.equals("S")) return "AGU";
        if(aminoAcid.equals("R")) return "AGA";
        if(aminoAcid.equals("W")) return "UGG";
        if(aminoAcid.equals("R")) return "AGG";
        if(aminoAcid.equals("L")) return "CUA";
        if(aminoAcid.equals("V")) return "GUA";
        if(aminoAcid.equals("L")) return "CUC";
        if(aminoAcid.equals("V")) return "GUC";
        if(aminoAcid.equals("L")) return "CUG";
        if(aminoAcid.equals("V")) return "GUG";
        if(aminoAcid.equals("L")) return "CUU";
        if(aminoAcid.equals("V")) return "GUU";
        if(aminoAcid.equals("P")) return "CCA";
        if(aminoAcid.equals("A")) return "GCA";
        if(aminoAcid.equals("P")) return "CCC";
        if(aminoAcid.equals("A")) return "GCC";
        if(aminoAcid.equals("P")) return "CCG";
        if(aminoAcid.equals("A")) return "GCG";
        if(aminoAcid.equals("P")) return "CCU";
        if(aminoAcid.equals("A")) return "GCU";
        if(aminoAcid.equals("H")) return "CAC";
        if(aminoAcid.equals("D")) return "GAC";
        if(aminoAcid.equals("H")) return "CAU";
        if(aminoAcid.equals("D")) return "GAU";
        if(aminoAcid.equals("Q")) return "CAA";
        if(aminoAcid.equals("E")) return "GAA";
        if(aminoAcid.equals("Q")) return "CAG";
        if(aminoAcid.equals("E")) return "GAG";
        if(aminoAcid.equals("R")) return "CGA";
        if(aminoAcid.equals("G")) return "GGA";
        if(aminoAcid.equals("R")) return "CGC";
        if(aminoAcid.equals("G")) return "GGC";
        if(aminoAcid.equals("R")) return "CGG";
        if(aminoAcid.equals("G")) return "GGG";
        if(aminoAcid.equals("R")) return "CGU";
        if(aminoAcid.equals("G")) return "GGU";

        return UNKNOWN;
    }

}
