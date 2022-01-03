package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.AminoAcids.AMINO_ACID_TO_CODON_MAP;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseStrandBases;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;

import java.util.List;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;

public final class ImpactTestUtils
{
    public static String generateAlt(final String ref)
    {
        String alt = "";
        for(int i = 0; i < ref.length(); ++i)
        {
            alt += getNextBase(ref.charAt(i));
        }

        return alt;
    }

    // basic short transcripts for testing
    public static MockRefGenome createMockGenome() { return createMockGenome(150); }

    public static MockRefGenome createMockGenome(int requiredBases)
    {
        final MockRefGenome refGenome = new MockRefGenome();
        String chr1Bases = generateRandomBases(requiredBases);
        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);
        return refGenome;
    }

    public static TranscriptData createPosTranscript()
    {
        // codons: 15-17, 18-20, 30-32, 33-35, 36-38, 39-50, 51-53 etc
        int[] exonStarts = { 10, 30, 50, 70, 90 };

        Integer codingStart = Integer.valueOf(15);
        Integer codingEnd = Integer.valueOf(75);

        return createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");
    }

    public static TranscriptData createNegTranscript()
    {
        // codons: 95-93, 92-90, 80-78, 77-75, 74-72, 71-60, 59-57 etc
        int[] exonStarts = { 10, 30, 50, 70, 90, 110, 130 };

        Integer codingStart = Integer.valueOf(15);
        Integer codingEnd = Integer.valueOf(95);

        return createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");
    }

    public static VariantData createSnv(int position, final String refBases)
    {
        String ref = refBases.substring(position, position + 1);
        String alt = getNextBase(ref);
        return new VariantData(CHR_1, position, ref, alt);
    }

    public static String getAminoAcidsCodons(final String aminoAcids, boolean reverseStrand)
    {
        String codonBases = "";

        // A E F on reverse becomes rev(F) rev(E) rev(A)
        for(int i = 0; i < aminoAcids.length(); ++i)
        {
            if(reverseStrand)
                codonBases = reverseStrandBases(getAminoAcidCodon(aminoAcids.charAt(i))) + codonBases;
            else
                codonBases += getAminoAcidCodon(aminoAcids.charAt(i));
        }

        return codonBases;
    }

    public static String getAminoAcidCodon(final char aminoAcid) { return getAminoAcidCodon(aminoAcid, 0); }

    public static String getAminoAcidCodon(final char aminoAcid, int index)
    {
        List<String> codons = AMINO_ACID_TO_CODON_MAP.get(String.valueOf(aminoAcid));

        if(codons == null || index >= codons.size())
            return "err";

        return codons.get(index);
    }
}
