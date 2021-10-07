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
import static com.hartwig.hmftools.pave.ImpactTestUtils.createSnv;

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
