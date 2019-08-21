package com.hartwig.hmftools.linx.neoepitope;

public class AminoAcidConverter
{
    public static String convertCodonToAminoAcid(final String codon)
    {
        if(codon.equals("UAA") || codon.equals("UAG") || codon.equals("UGA")) return "_";

        if(codon.equals("UCA")) return "S";
        if(codon.equals("AUA")) return "I";
        if(codon.equals("UCC")) return "S";
        if(codon.equals("AUC")) return "I";
        if(codon.equals("UCG")) return "S";
        if(codon.equals("AUU")) return "I";
        if(codon.equals("UCU")) return "S";
        if(codon.equals("AUG")) return "M";
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

        return "";
    }
}
