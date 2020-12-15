package com.hartwig.hmftools.common.codon;

public class Codons {

    public static String asCodonString(String dna) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < dna.length() - 2; i += 3) {
            builder.append(asCondon(dna.substring(i, i + 3)));
        }

        return builder.toString();
    }

    public static char asCondon(String dna) {
        assert (dna.length() == 3);
        switch (dna) {
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
            case "ATG":
                return 'M';
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
            case "TAA":
            case "TAG":
                return 'X';
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
            case "TGA":
                return 'X';
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

        return '.';
    }

}
