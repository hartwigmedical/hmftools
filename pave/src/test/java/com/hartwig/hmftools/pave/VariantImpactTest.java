package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Codons.START_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.STOP_CODON_1;
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
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FIVE_PRIME_UTR;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INTRONIC;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.NON_CODING_TRANSCRIPT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.START_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.THREE_PRIME_UTR;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.UPSTREAM_GENE;
import static com.hartwig.hmftools.pave.ImpactClassifier.checkStopStartCodons;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.codon.AminoAcids;
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

        int[] exonStarts = { 100, 200, 300, 400 };
        Integer codingStart = new Integer(125);
        Integer codingEnd = new Integer(425);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 50, codingStart, codingEnd, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // pre-gene
        int pos = 50;
        VariantData var = createSnv(pos, refBases);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(UPSTREAM_GENE, impact.topEffect());

        // 5' UTR
        pos = 120;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(FIVE_PRIME_UTR, impact.topEffect());

        // 3' UTR
        pos = 440;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(THREE_PRIME_UTR, impact.topEffect());

        // intronic
        pos = 175;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(INTRONIC, impact.topEffect());

        TranscriptData transDataNegStrand = createTransExons(
                GENE_ID_2, TRANS_ID_2, NEG_STRAND, exonStarts, 50, codingStart, codingEnd, false, "");

        // pre-gene
        pos = 490;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataNegStrand);
        assertEquals(UPSTREAM_GENE, impact.topEffect());

        // intronic
        pos = 375;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataNegStrand);
        assertEquals(INTRONIC, impact.topEffect());

        // non-coding exonic
        TranscriptData transDataNonCoding = createTransExons(
                GENE_ID_2, TRANS_ID_2, NEG_STRAND, exonStarts, 50, null, null, false, "");

        // intronic
        pos = 160;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataNonCoding);
        assertEquals(INTRONIC, impact.topEffect());

        // exonic has special classification
        pos = 220;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataNonCoding);
        assertEquals(NON_CODING_TRANSCRIPT, impact.topEffect());
    }

    private VariantData createSnv(int position, final String refBases)
    {
        String ref = refBases.substring(position, position + 1);
        String alt = getNextBase(ref);
        return new VariantData(CHR_1, position, ref, alt);
    }

    @Test
    public void testSynonymousMissenseImpacts()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        int[] exonStarts = { 0, 100, 200 };

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
        assertEquals(SYNONYMOUS, impact.topEffect());

        // now missense
        pos = 41;
        ref = chr1Bases.substring(pos, pos + 1);
        alt = "T";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(MISSENSE, impact.topEffect());

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
        assertEquals(SYNONYMOUS, impact.topEffect());

        // test for an MNV spanning 3 codons - first in the last codon pos at 42 then all the next and 2 into the final one at 46
        pos = 42;
        ref = chr1Bases.substring(pos, pos + 5);
        alt = "G" + "ATT" + "C";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(SYNONYMOUS, impact.topEffect());

        // now missense by change middle codon
        alt = "G" + "AAT" + "C";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(MISSENSE, impact.topEffect());

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
        assertEquals(SYNONYMOUS, impact.topEffect());

        // now missense
        pos = 243;
        ref = chr1Bases.substring(pos, pos + 1);
        alt = "T";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNegStrand);
        assertEquals(MISSENSE, impact.topEffect());

        // test a MNV spanning 3 codons as before
        pos = 241;
        ref = chr1Bases.substring(pos, pos + 5);
        alt = reverseStrandBases("G" + "ATT" + "C");
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNegStrand);
        assertEquals(SYNONYMOUS, impact.topEffect());

        // and missense
        alt = reverseStrandBases("G" + "AAT" + "C");
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNegStrand);
        assertEquals(MISSENSE, impact.topEffect());
    }

    @Test
    public void testInframeIndelImpacts()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(300);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = { 0, 100, 200 };
        Integer codingStart = new Integer(10);
        Integer codingEnd = new Integer(250);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 80, codingStart, codingEnd, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // inframe conservative deletion of 2 codons (ie 3 -> 1)
        int pos = 48; // 3rd base of codon
        String ref = chr1Bases.substring(pos, pos + 7);
        String alt = chr1Bases.substring(pos, pos + 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(INFRAME_DELETION, impact.topEffect());

        // inframe disruptive (spans first and last codon boundaries) deletion of 3 codons (ie 4 -> 1)
        pos = 47; // first base of codon
        ref = chr1Bases.substring(pos, pos + 10);
        alt = chr1Bases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(INFRAME_DELETION, impact.topEffect());

        // conservative insertion of 2 codons
        pos = 46; // first base of codon
        ref = chr1Bases.substring(pos, pos + 1);
        alt = ref + "AGAGAG";
        var = new VariantData(CHR_1, pos, ref, alt);

        // TODO

        impact = classifier.classifyVariant(var, transDataPosStrand);
        // assertEquals(INFRAME_INSERTION, impact.consequence());

        // disruptive insertion of 3 codons
        pos = 45; // first base of codon
        ref = chr1Bases.substring(pos, pos + 1);
        alt = ref + generateRandomBases(9);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        // assertEquals(INFRAME_INSERTION, impact.consequence());
    }

    @Test
    public void testNonsenseImpacts()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(300);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = { 0, 100, 200 };
        Integer codingStart = new Integer(10);
        Integer codingEnd = new Integer(250);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 80, codingStart, codingEnd, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // stop gained by MNV
        int pos = 46; // first base of codon
        String ref = chr1Bases.substring(pos, pos + 3);
        String alt = STOP_CODON_1;
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(STOP_GAINED, impact.topEffect());
    }

    @Test
    public void testStopStartLostGained()
    {
        // stop gained
        String refAminoAcids = "SILFT";
        String altAminoAcids = "SI" + STOP_AMINO_ACID + "FT";

        assertEquals(STOP_GAINED, checkStopStartCodons(10, refAminoAcids, altAminoAcids));

        // stop lost
        refAminoAcids = "SILF" + STOP_AMINO_ACID;
        altAminoAcids = "SILFTPW";

        assertEquals(STOP_LOST, checkStopStartCodons(10, refAminoAcids, altAminoAcids));

        // start lost
        refAminoAcids = START_AMINO_ACID + "SILFT";
        altAminoAcids = "WSILFT";

        assertEquals(START_LOST, checkStopStartCodons(1, refAminoAcids, altAminoAcids));

        assertEquals(null, checkStopStartCodons(2, refAminoAcids, altAminoAcids));

    }
}
