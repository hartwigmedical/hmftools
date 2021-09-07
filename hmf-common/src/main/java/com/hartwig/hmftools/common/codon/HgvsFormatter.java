package com.hartwig.hmftools.common.codon;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

// TODO - reconcile with AminoAcidFunctions
// Do note that according to http://www.hgvs.org/mutnomen/standards.html there are also 2 "dubious" mappings (not included below).
public final class HgvsFormatter
{

    private static final String HGVS_CODING_PREFIX_TO_REMOVE = "c.";
    private static final String HGVS_PROTEIN_PREFIX_TO_REMOVE = "p.";

    private static final Map<String, String> ONE_TO_THREE_AMINO_ACID_MAPPING = Maps.newHashMap();

    static
    {
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("A", "Ala"); // Alanine (GCA, GCC, GCG, GCT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("C", "Cys"); // Cysteine (TGC, TGT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("D", "Asp"); // Aspartic acid (GAC, GAT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("E", "Glu"); // Glutamic acid (GAA, GAG)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("F", "Phe"); // Phenylalanine (TTC, TTT)

        ONE_TO_THREE_AMINO_ACID_MAPPING.put("G", "Gly"); // Glycine (GGA, GGC, GGG, GGT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("H", "His"); // Histidine (CAC, CAT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("I", "Ile"); // Isoleucine (ATA, ATC, ATT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("K", "Lys"); // Lysine (AAA, AAG)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("L", "Leu"); // Leucine (TTA, TTG, CTA, CTC, CTG, CTT)

        ONE_TO_THREE_AMINO_ACID_MAPPING.put("M", "Met"); // Methionine (ATG)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("N", "Asn"); // Asparagine (AAC, AAT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("P", "Pro"); // Proline (CCA, CCC, CCG, CCT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("Q", "Gln"); // Glutamine (CAA, CAG)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("R", "Arg"); // Arginine (AGA, AGG, CGA, CGC, CGG, CGT)

        ONE_TO_THREE_AMINO_ACID_MAPPING.put("S", "Ser"); // Serine (AGC, AGT, TCA, TCC, TCG, TCT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("T", "Thr"); // Threonine (ACA, ACC, ACG, ACT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("V", "Val"); // Valine (GTA, GTC, GTG, GTT)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("W", "Trp"); // Tryptophan (TGG)
        ONE_TO_THREE_AMINO_ACID_MAPPING.put("Y", "Tyr"); // Tyrosine (TAC, TAT)
    }

    @NotNull
    public static String formattedHgvsCoding(final String hgvsCoding)
    {
        if(hgvsCoding.startsWith(HGVS_CODING_PREFIX_TO_REMOVE))
            return hgvsCoding.substring(HGVS_CODING_PREFIX_TO_REMOVE.length());
        else
            return hgvsCoding;
    }

    @NotNull
    public static String formattedHgvsProtein(final String hgvsProteinRaw)
    {
        String hgvsProtein = hgvsProteinRaw;

        if(hgvsProtein.startsWith(HGVS_PROTEIN_PREFIX_TO_REMOVE))
            hgvsProtein = hgvsProtein.substring(HGVS_PROTEIN_PREFIX_TO_REMOVE.length());

        for(Map.Entry<String, String> mappingEntry : ONE_TO_THREE_AMINO_ACID_MAPPING.entrySet())
        {
            hgvsProtein = hgvsProtein.replaceAll(mappingEntry.getValue(), mappingEntry.getKey());
        }

        return hgvsProtein;
    }
}
