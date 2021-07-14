package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseStrandBases;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INTRON_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.MISSENSE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SYNONYMOUS_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.UPSTREAM_GENE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.UTR_VARIANT;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.Set;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class VariantImpactTest
{

    @Test
    public void testNonCodingImpacts()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String refBases = generateRandomBases(500);

        refGenome.RefGenomeMap.put(CHR_1, refBases);

        int[] exonStarts = {100, 200, 300, 400};
        Integer codingStart = new Integer(125);
        Integer codingEnd = new Integer(425);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 50, codingStart, codingEnd, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // pre-gene
        int pos = 50;
        VariantData var = createSnv(pos, refBases);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(UPSTREAM_GENE_VARIANT, impact.Consequence);

        // 5' UTR
        pos = 120;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(UTR_VARIANT, impact.Consequence);

        // 3' UTR
        pos = 440;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(UTR_VARIANT, impact.Consequence);

        // intronic
        pos = 175;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(INTRON_VARIANT, impact.Consequence);

        TranscriptData transDataNegStrand = createTransExons(
                GENE_ID_2, TRANS_ID_2, NEG_STRAND, exonStarts, 50, codingStart, codingEnd, false, "");

        // pre-gene
        pos = 490;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataNegStrand);
        assertEquals(UPSTREAM_GENE_VARIANT, impact.Consequence);

        // intronic
        pos = 375;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataNegStrand);
        assertEquals(INTRON_VARIANT, impact.Consequence);
    }

    private VariantData createSnv(int position, final String refBases)
    {
        String ref = refBases.substring(position, position + 1);
        String alt = getNextBase(ref);
        return new VariantData(CHR_1, position, ref, alt);
    }

    @Test
    public void testSynonymousImpacts()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        int[] exonStarts = {0, 100, 200};

        // codons start on at 10, 13, 16 etc
        Integer codingStart = new Integer(10);
        Integer codingEnd = new Integer(250);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 80, codingStart, codingEnd, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        String chr1Bases = generateRandomBases(300);

        // set the specific AAs for a region of this mock ref genome
        // S: TCA -> TCG - so last base of codon can change
        // I: ATC -> ATT
        // L: TTA -> CTA - first base changes

        //                     40      43      46
        String aminoAcidSeq = "TCA" + "ATC" + "TTA";
        chr1Bases = chr1Bases.substring(0, 40) + aminoAcidSeq + chr1Bases.substring(49);
        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        // test SNVs at the 1st and 3rd codon bases

        // last base of a codon changes
        int codonPos = 40;
        String codon = chr1Bases.substring(codonPos, codonPos + 3);
        String aminoAcid = AminoAcids.findAminoAcidForCodon(codon);

        String alt = "G";
        String synCodon = chr1Bases.substring(codonPos, codonPos + 2) + alt;
        assertTrue(aminoAcid.equals(AminoAcids.findAminoAcidForCodon(synCodon)));

        int pos = 42;
        String ref = chr1Bases.substring(pos, pos + 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(SYNONYMOUS_VARIANT, impact.Consequence);

        // first base of a codon changes
        codonPos = 46;
        codon = chr1Bases.substring(codonPos, codonPos + 3);
        aminoAcid = AminoAcids.findAminoAcidForCodon(codon);

        alt = "C";
        synCodon = alt + chr1Bases.substring(codonPos + 1, codonPos + 3);
        assertTrue(aminoAcid.equals(AminoAcids.findAminoAcidForCodon(synCodon)));

        pos = 46;
        ref = chr1Bases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(SYNONYMOUS_VARIANT, impact.Consequence);


        // test for an MNV spanning 3 codons - first in the last codon pos at 42 then all the next and 2 into the final one at 46
        pos = 42;
        ref = chr1Bases.substring(pos, pos + 5);
        alt = "G" + "ATT" + "C";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(SYNONYMOUS_VARIANT, impact.Consequence);

        // test reverse strand
        TranscriptData transDataNegStrand = createTransExons(
                GENE_ID_2, TRANS_ID_2, NEG_STRAND, exonStarts, 80, codingStart, codingEnd, false, "");

        // coding starts at 250 so codons start at 250, 247, 244 etc
        // still change the AA sequence S, I then L - for the range 239-241, 242-244 and 245-247
        String aminoAcidSeqRev = reverseStrandBases(aminoAcidSeq);
        chr1Bases = chr1Bases.substring(0, 239) + aminoAcidSeqRev + chr1Bases.substring(248);
        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        codonPos = 242;
        codon = chr1Bases.substring(codonPos, codonPos + 3);
        aminoAcid = AminoAcids.findAminoAcidForCodon(reverseStrandBases(codon));

        // change last base of a codon, which is the first base
        alt = "A";
        synCodon = alt + chr1Bases.substring(codonPos + 1, codonPos + 3);
        assertTrue(aminoAcid.equals(AminoAcids.findAminoAcidForCodon(reverseStrandBases(synCodon))));

        pos = 242;
        ref = chr1Bases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNegStrand);
        assertEquals(SYNONYMOUS_VARIANT, impact.Consequence);

        // test a MNV spanning 3 codons as before
        pos = 241;
        ref = chr1Bases.substring(pos, pos + 5);
        alt = reverseStrandBases("G" + "ATT" + "C");
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNegStrand);
        assertEquals(SYNONYMOUS_VARIANT, impact.Consequence);
    }

    /*
    @Test
    public void testMissenseImpacts()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(300);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = {0, 100, 200};
        Integer codingStart = new Integer(10);
        Integer codingEnd = new Integer(250);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 80, codingStart, codingEnd, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        int pos = 28;
        String ref = chr1Bases.substring(pos, pos + 1);
        String alt = "G";
        VariantData var = new VariantData(CHR_1, 28, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(MISSENSE_VARIANT, impact.Consequence);

    }
     */

    private String findMissenseBase()
    {
        return "";
    }
}
