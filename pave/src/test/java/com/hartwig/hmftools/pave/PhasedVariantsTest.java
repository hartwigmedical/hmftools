package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.PHASED_SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createMockGenome;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateAlt;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateTestBases;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import org.junit.Test;

public class PhasedVariantsTest
{
    /* scenarios:
        - DEL and DEL causing inframe
        - DEL and DEL remaining frameshift
        - DEL and INS causing missense
        - DEL and INS causing synonymous
        - 3x DEL/INS causing inframe insert or delete or missense
        - SNVs before and after a set of phaseable INDELs
    */

    private MockRefGenome mRefGenome = createMockGenome(201);
    private String mRefBases = mRefGenome.RefGenomeMap.get(CHR_1);
    private ImpactClassifier mClassifier = new ImpactClassifier(mRefGenome);

    @Test
    public void testInframeDelPair()
    {
        TranscriptData transDataPos = createTransExons(
            GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] {0, 100}, 80, 10, 190, false, "");

        // 2 DELs making inframe DEL
        int pos = 20;

        String ref = mRefBases.substring(pos, pos + 3);
        String alt = mRefBases.substring(pos, pos + 1);
        VariantData var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact1 = classifyVariant(var1, transDataPos);

        assertTrue(impact1.codingContext().IsFrameShift);
        assertTrue(impact1.hasEffect(FRAMESHIFT));

        pos = 25;
        ref = mRefBases.substring(pos, pos + 2);
        alt = mRefBases.substring(pos, pos + 1);
        VariantData var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact2 = classifyVariant(var2, transDataPos);

        assertTrue(impact2.codingContext().IsFrameShift);
        assertTrue(impact2.hasEffect(FRAMESHIFT));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertTrue(impact1.phasedFrameshift());
        assertFalse(impact1.codingContext().IsFrameShift);
        assertFalse(impact1.hasEffect(FRAMESHIFT));
        assertTrue(impact1.hasEffect(PHASED_INFRAME_DELETION));
        assertFalse(impact2.codingContext().IsFrameShift);
        assertFalse(impact2.hasEffect(FRAMESHIFT));
        assertTrue(impact2.hasEffect(PHASED_INFRAME_DELETION));

        // as before but with an extra base so remaining out of frame
        var1.getImpacts().clear();
        impact1 = classifyVariant(var1, transDataPos);

        assertTrue(impact1.codingContext().IsFrameShift);
        assertTrue(impact1.hasEffect(FRAMESHIFT));

        pos = 25;
        ref = mRefBases.substring(pos, pos + 3);
        alt = mRefBases.substring(pos, pos + 1);
        var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        impact2 = classifyVariant(var2, transDataPos);

        assertTrue(impact2.codingContext().IsFrameShift);
        assertTrue(impact2.hasEffect(FRAMESHIFT));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertTrue(impact1.codingContext().IsFrameShift);
        assertTrue(impact1.hasEffect(FRAMESHIFT));

        assertTrue(impact2.codingContext().IsFrameShift);
        assertTrue(impact2.hasEffect(FRAMESHIFT));

        // 2 DELs which overlap ref codon bases but can form a frameshift
        // codons 13-15, 16-18, 19-21, 22-24 etc
        pos = 16;
        ref = mRefBases.substring(pos, pos + 3);
        alt = mRefBases.substring(pos, pos + 1); // deletes 17 & 18;
        var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        impact1 = classifyVariant(var1, transDataPos);
        assertEquals(16, impact1.proteinContext().refCodingBaseStart());
        assertEquals(21, impact1.proteinContext().refCodingBaseEnd());

        assertTrue(impact1.codingContext().IsFrameShift);
        assertTrue(impact1.hasEffect(FRAMESHIFT));

        pos = 20;
        ref = mRefBases.substring(pos, pos + 2);
        alt = mRefBases.substring(pos, pos + 1); // deletes 21 being the 3rd codon base which they share
        var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        impact2 = classifyVariant(var2, transDataPos);
        assertEquals(19, impact2.proteinContext().refCodingBaseStart());
        assertEquals(24, impact2.proteinContext().refCodingBaseEnd());

        assertTrue(impact2.codingContext().IsFrameShift);
        assertTrue(impact2.hasEffect(FRAMESHIFT));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertTrue(impact1.phasedFrameshift());
        assertFalse(impact1.codingContext().IsFrameShift);
        assertFalse(impact1.hasEffect(FRAMESHIFT));
        assertTrue(impact1.hasEffect(PHASED_INFRAME_DELETION));
        assertFalse(impact2.codingContext().IsFrameShift);
        assertFalse(impact2.hasEffect(FRAMESHIFT));
        assertTrue(impact2.hasEffect(PHASED_INFRAME_DELETION));
    }

    @Test
    public void testInframeDelInsPair()
    {
        TranscriptData transDataNeg = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, new int[] { 0, 100 }, 80, 10, 190, false, "");

        // DEL and INS making inframe insert
        int pos = 20;

        String ref = mRefBases.substring(pos, pos + 3);
        String alt = mRefBases.substring(pos, pos + 1);
        VariantData var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact1 = classifyVariant(var1, transDataNeg);

        assertTrue(impact1.codingContext().IsFrameShift);
        assertTrue(impact1.hasEffect(FRAMESHIFT));

        pos = 25;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "AAAAAAAA";
        VariantData var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact2 = classifyVariant(var2, transDataNeg);

        assertTrue(impact2.codingContext().IsFrameShift);
        assertTrue(impact2.hasEffect(FRAMESHIFT));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertTrue(impact1.phasedFrameshift());
        assertFalse(impact1.codingContext().IsFrameShift);
        assertFalse(impact1.hasEffect(FRAMESHIFT));
        assertTrue(impact1.hasEffect(PHASED_INFRAME_INSERTION));
        assertFalse(impact2.codingContext().IsFrameShift);
        assertFalse(impact2.hasEffect(FRAMESHIFT));
        assertTrue(impact2.hasEffect(PHASED_INFRAME_INSERTION));

        // again but remaining out-of-frame
        var1.getImpacts().clear();

        impact1 = classifyVariant(var1, transDataNeg);

        pos = 25;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "AAAAAAA";
        var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        impact2 = classifyVariant(var2, transDataNeg);

        assertTrue(impact2.codingContext().IsFrameShift);
        assertTrue(impact2.hasEffect(FRAMESHIFT));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertTrue(impact1.codingContext().IsFrameShift);
        assertTrue(impact1.hasEffect(FRAMESHIFT));
        assertTrue(impact2.codingContext().IsFrameShift);
        assertTrue(impact2.hasEffect(FRAMESHIFT));
    }

    @Test
    public void testDelSnvOverlapPair()
    {
        TranscriptData transDataPos = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] { 0, 100 }, 80, 10, 190, false, "");

        // DEL and SNV, overlapping in the same codon (ala EGFR)
        int pos = 26;

        String ref = mRefBases.substring(pos, pos + 10);
        String alt = mRefBases.substring(pos, pos + 1);

        VariantData var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact1 = classifyVariant(var1, transDataPos);

        assertTrue(impact1.hasEffect(INFRAME_DELETION));
        assertEquals("IDRS", impact1.proteinContext().RefAminoAcids);
        assertEquals("M", impact1.proteinContext().AltAminoAcids);
        pos = 36;
        ref = mRefBases.substring(pos, pos + 1);
        alt =  "T";

        VariantData var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact2 = classifyVariant(var2, transDataPos);

        assertTrue(impact2.hasEffect(SYNONYMOUS));
        assertEquals("S", impact2.proteinContext().AltAminoAcids);

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertTrue(impact1.phasedFrameshift());
        assertFalse(impact1.hasEffect(FRAMESHIFT));
        assertTrue(impact1.hasEffect(PHASED_INFRAME_DELETION));
        assertFalse(impact2.codingContext().IsFrameShift);
        assertFalse(impact2.hasEffect(FRAMESHIFT));
        assertTrue(impact2.hasEffect(PHASED_INFRAME_DELETION));

        // again but with the DEL conservative and not overlapping with the SNV despite its ref codon bases overlapping
        // 13-18 are deleted, SNV is at 20
        pos = 12;
        ref = mRefBases.substring(pos, pos + 7);
        alt = mRefBases.substring(pos, pos + 1);

        var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        impact1 = classifyVariant(var1, transDataPos);

        assertTrue(impact1.hasEffect(INFRAME_DELETION));
        assertEquals("SIDR", impact1.proteinContext().RefAminoAcids);
        assertEquals("SR", impact1.proteinContext().AltAminoAcids);

        pos = 20;
        ref = mRefBases.substring(pos, pos + 1);
        alt =  "T";

        var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        impact2 = classifyVariant(var2, transDataPos);

        assertTrue(impact2.hasEffect(MISSENSE));
        assertEquals("R", impact2.proteinContext().RefAminoAcids);
        assertEquals("L", impact2.proteinContext().AltAminoAcids);

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertFalse(impact1.phasedFrameshift());
    }

    @Test
    public void testMissenseSynonymousDelInsPair()
    {
        TranscriptData transDataNeg = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] { 0, 100 }, 80, 10, 190, false, "");

        // DEL and INS making missense and symonymous effects
        int pos = 20;

        String ref = mRefBases.substring(pos, pos + 1);
        String alt = ref + "AGCT";
        VariantData var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact1 = classifyVariant(var1, transDataNeg);

        assertTrue(impact1.codingContext().IsFrameShift);
        assertTrue(impact1.hasEffect(FRAMESHIFT));

        pos = 25;
        ref = mRefBases.substring(pos, pos + 5);
        alt = mRefBases.substring(pos, pos + 1);
        VariantData var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact2 = classifyVariant(var2, transDataNeg);

        assertTrue(impact2.codingContext().IsFrameShift);
        assertTrue(impact2.hasEffect(FRAMESHIFT));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertTrue(impact1.phasedFrameshift());
        assertFalse(impact1.codingContext().IsFrameShift);
        assertFalse(impact1.hasEffect(FRAMESHIFT));
        assertTrue(impact1.hasEffect(PHASED_MISSENSE));
        assertFalse(impact2.codingContext().IsFrameShift);
        assertFalse(impact2.hasEffect(FRAMESHIFT));
        assertTrue(impact2.hasEffect(PHASED_MISSENSE));
        assertEquals("p.Ser5_Ile6delinsAlaIle", impact1.proteinContext().Hgvs);

        // now synonynous
        pos = 20;
        ref = mRefBases.substring(pos, pos + 5);
        alt = mRefBases.substring(pos, pos + 1);
        var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(2, "", "", 0);
        impact1 = classifyVariant(var1, transDataNeg);

        assertTrue(impact1.codingContext().IsFrameShift);
        assertTrue(impact1.hasEffect(FRAMESHIFT));

        pos = 28;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "ATCG";
        var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(2, "", "", 0);
        impact2 = classifyVariant(var2, transDataNeg);

        assertTrue(impact2.codingContext().IsFrameShift);
        assertTrue(impact2.hasEffect(FRAMESHIFT));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertTrue(impact1.phasedFrameshift());
        assertFalse(impact1.codingContext().IsFrameShift);
        assertFalse(impact1.hasEffect(FRAMESHIFT));
        assertTrue(impact1.hasEffect(PHASED_SYNONYMOUS));
        assertFalse(impact2.codingContext().IsFrameShift);
        assertFalse(impact2.hasEffect(FRAMESHIFT));
        assertTrue(impact2.hasEffect(PHASED_SYNONYMOUS));
        assertEquals("p.Arg4_Arg8=", impact1.proteinContext().Hgvs);
    }

    @Test
    public void testMultiIndel()
    {
        TranscriptData transDataNeg = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, new int[] { 0, 100 }, 80, 10, 190, false, "");

        // DEL and INS making missense and symonymous effects
        // SNVs at either end are ignored
        int pos = 14;
        String ref = mRefBases.substring(pos, pos + 1);
        String alt = generateAlt(ref);
        VariantData varIgnore1 = new VariantData(CHR_1, pos, ref, alt);
        varIgnore1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impactIgnore1 = classifyVariant(varIgnore1, transDataNeg);
        assertTrue(impactIgnore1.hasEffect(MISSENSE));

        pos = 20;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "AG";
        VariantData var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact1 = classifyVariant(var1, transDataNeg);

        assertTrue(impact1.codingContext().IsFrameShift);
        assertTrue(impact1.hasEffect(FRAMESHIFT));

        pos = 25;
        ref = mRefBases.substring(pos, pos + 8);
        alt = mRefBases.substring(pos, pos + 1);
        VariantData var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact2 = classifyVariant(var2, transDataNeg);

        assertTrue(impact2.codingContext().IsFrameShift);
        assertTrue(impact2.hasEffect(FRAMESHIFT));

        pos = 35;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "GG";
        VariantData var3 = new VariantData(CHR_1, pos, ref, alt);
        var3.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact3 = classifyVariant(var3, transDataNeg);

        assertTrue(impact3.codingContext().IsFrameShift);
        assertTrue(impact3.hasEffect(FRAMESHIFT));

        pos = 40;
        ref = mRefBases.substring(pos, pos + 1);
        alt = generateAlt(ref);
        VariantData varIgnore2 = new VariantData(CHR_1, pos, ref, alt);
        varIgnore2.setVariantDetails(1, "", "", 0);
        VariantTransImpact impactIgnore2 = classifyVariant(varIgnore2, transDataNeg);
        assertTrue(impactIgnore2.hasEffect(SYNONYMOUS));

        pos = 46;
        ref = mRefBases.substring(pos, pos + 1);
        alt = generateAlt(ref);
        VariantData varIgnore3 = new VariantData(CHR_1, pos, ref, alt);
        varIgnore3.setVariantDetails(1, "", "", 0);
        VariantTransImpact impactIgnore3 = classifyVariant(varIgnore3, transDataNeg);
        assertTrue(impactIgnore3.hasEffect(SYNONYMOUS));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertFalse(impact1.codingContext().IsFrameShift);
        assertFalse(impact1.hasEffect(FRAMESHIFT));
        assertTrue(impact1.phasedFrameshift());
        assertTrue(impact1.hasEffect(PHASED_INFRAME_DELETION));
        assertTrue(impact2.hasEffect(PHASED_INFRAME_DELETION));
        assertTrue(impact3.hasEffect(PHASED_INFRAME_DELETION));

        assertFalse(impactIgnore1.phasedFrameshift());
        assertFalse(impactIgnore2.phasedFrameshift());
        assertFalse(impactIgnore3.phasedFrameshift());
        assertFalse(impactIgnore1.hasEffect(PHASED_INFRAME_DELETION));
        assertFalse(impactIgnore2.hasEffect(PHASED_INFRAME_DELETION));
        assertFalse(impactIgnore3.hasEffect(PHASED_INFRAME_DELETION));
    }

    @Test
    public void testRefCondonOverlap()
    {
        // variants remain as-is since do not overlap each other's impacted ref or alt codons

        TranscriptData transData = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] { 0, 100 }, 80, 10, 190, false, "");

        // positions:      012345678901234567890123     456789012345     67890
        String refBases = "XCCCCCCCCCCAAAAAAAAAACCC" + "GGTGGCGGCGGC" + "AAAAATTTTTAAAAATTTTT";

        int pos = 26;
        String ref = refBases.substring(pos, pos + 7);
        String alt = refBases.substring(pos, pos + 1);
        VariantData varIgnore1 = new VariantData(CHR_1, pos, ref, alt);
        varIgnore1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact1 = classifyVariant(varIgnore1, transData);
        assertTrue(impact1.hasEffect(INFRAME_DELETION));

        pos = 35;
        ref = refBases.substring(pos, pos + 1);
        alt = "T";
        VariantData var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact2 = classifyVariant(var1, transData);
        assertTrue(impact2.hasEffect(MISSENSE));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertFalse(impact1.phasedFrameshift());
        assertTrue(impact1.hasEffect(INFRAME_DELETION));
        assertFalse(impact2.phasedFrameshift());
        assertTrue(impact2.hasEffect(MISSENSE));
    }

    @Test
    public void testSpliceEdgeOverlap()
    {
        TranscriptData transData = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] { 0, 20 }, 12, 1, 25, false, "");

        //                              intron
        //                           10        20        30        40
        // positions:      01234567890123456789012345678901234567890
        String refBases = "XATGCATGGGCAGTTTTTTTTTTTTTAAAAATTTTTAAAAATTTTT";

        mRefGenome.RefGenomeMap.put(CHR_2, refBases);

        // start
        int pos = 8;
        String ref = refBases.substring(pos, pos + 1);
        String alt = ref + "A";
        VariantData var1 = new VariantData(CHR_2, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact1 = classifyVariant(var1, transData);
        assertTrue(impact1.hasEffect(FRAMESHIFT));

        pos = 9;
        ref = refBases.substring(pos, pos + 1);
        alt = "T";
        VariantData var2 = new VariantData(CHR_2, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact2 = classifyVariant(var2, transData);
        assertTrue(impact2.hasEffect(SYNONYMOUS));

        pos = 11;
        ref = refBases.substring(pos, pos + 9); // AGGTAACTT>A
        alt = refBases.substring(pos, pos + 1);
        VariantData var3 = new VariantData(CHR_2, pos, ref, alt);
        var3.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact3 = classifyVariant(var3, transData);
        assertTrue(impact3.hasEffect(FRAMESHIFT));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertEquals("p.Gln4Ser", impact1.hgvsProtein());
        assertTrue(impact1.phasedFrameshift());
        assertTrue(impact1.hasEffect(PHASED_MISSENSE));
        assertTrue(impact2.hasEffect(PHASED_MISSENSE));
        assertTrue(impact3.hasEffect(PHASED_MISSENSE));
    }

    @Test
    public void testMixedOverlaps()
    {
        TranscriptData transData = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] { 0, 100 }, 80, 10, 190, false, "");
        // positions:      01234567890123456789012     34567890123456789012345678901234567
        String refBases = "XCCCCCCCCCCAAAAAAAAAACC" + "TTGAGAAGTCTGCCAGCACTGAGAGGAAAATTAAC" + "AAAAATTTTTAAAAATTTTT";

        int pos = 23;
        String ref = refBases.substring(pos, pos + 8);
        String alt = refBases.substring(pos, pos + 1);
        VariantData var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact1 = classifyVariant(var1, transData);
        assertTrue(impact1.hasEffect(FRAMESHIFT));

        pos = 34;
        ref = refBases.substring(pos, pos + 1);
        alt = "T";
        VariantData var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact2 = classifyVariant(var2, transData);
        assertTrue(impact2.hasEffect(SYNONYMOUS));

        pos = 35;
        ref = refBases.substring(pos, pos + 1);
        alt = ref + "T";
        VariantData var3 = new VariantData(CHR_1, pos, ref, alt);
        var3.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact3 = classifyVariant(var3, transData);
        assertTrue(impact3.hasEffect(FRAMESHIFT));

        pos = 36;
        ref = refBases.substring(pos, pos + 13);
        alt = refBases.substring(pos, pos + 1);
        VariantData var4 = new VariantData(CHR_1, pos, ref, alt);
        var4.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact4 = classifyVariant(var4, transData);
        assertTrue(impact4.hasEffect(INFRAME_DELETION));

        // the last variant doesn't overlap
        pos = 53;
        ref = refBases.substring(pos, pos + 1);
        alt = "A";
        VariantData var5 = new VariantData(CHR_1, pos, ref, alt);
        var5.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact5 = classifyVariant(var5, transData);
        assertTrue(impact5.hasEffect(SYNONYMOUS));

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertFalse(impact2.codingContext().IsFrameShift);
        assertFalse(impact2.hasEffect(FRAMESHIFT));
        assertTrue(impact2.phasedFrameshift());
        assertTrue(impact1.hasEffect(PHASED_INFRAME_DELETION));
        assertTrue(impact2.hasEffect(PHASED_INFRAME_DELETION));
        assertTrue(impact3.hasEffect(PHASED_INFRAME_DELETION));
        assertTrue(impact4.hasEffect(PHASED_INFRAME_DELETION));
        assertTrue(impact5.hasEffect(SYNONYMOUS));
    }

    @Test
    public void testMixedOverlaps2()
    {
        TranscriptData transData = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, new int[] { 0 }, 80, 5, 70, false, "");

        // positions:      0123456789012345678901234567890123456789012345678901234567
        //                 0         10        20        30        40        50
        String refBases = "XCCCCATGAAAAAAAAAACCCCCCCCCCGGGGGGGAAGCTGATGTTCAGGAGTG" + generateTestBases(40);
        mRefGenome.RefGenomeMap.put(CHR_1, refBases);

        // pos(12:50384540-50384542) variant(GA>G)"
        int pos = 40;
        String ref = refBases.substring(pos, pos + 2);
        String alt = ref.substring(0, 1);
        VariantData var1 = new VariantData(CHR_1, pos, ref, alt);
        var1.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact1 = classifyVariant(var1, transData);

        // pos(12:50384543) variant(G>C)
        pos = 43;
        ref = refBases.substring(pos, pos + 1);
        alt = "C";
        VariantData var2 = new VariantData(CHR_1, pos, ref, alt);
        var2.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact2 = classifyVariant(var2, transData);

        // pos(12:50384545-50384546) variant(T>TCTCAAGA)
        pos = 45;
        ref = refBases.substring(pos, pos + 1);
        alt = ref + "CTCAAGA";
        VariantData var3 = new VariantData(CHR_1, pos, ref, alt);
        var3.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact3 = classifyVariant(var3, transData);

        assertTrue(impact3.codingContext().IsFrameShift);
        assertTrue(impact3.hasEffect(FRAMESHIFT));

        // pos(12:50384548-50384549) variant(G>GCTT)
        pos = 48;
        ref = refBases.substring(pos, pos + 1);
        alt = ref + "CTT";
        VariantData var4 = new VariantData(CHR_1, pos, ref, alt);
        var4.setVariantDetails(1, "", "", 0);
        VariantTransImpact impact4 = classifyVariant(var4, transData);

        mClassifier.processPhasedVariants(NO_LOCAL_PHASE_SET);

        assertEquals("CTGATGTTCAGG", impact4.proteinContext().RefCodonBases);
        assertEquals("CTGTCTTCTCAAGACAGCTTG", impact4.proteinContext().AltCodonBases);
        assertEquals("p.Pro8_His10delinsGlnAlaValLeuArgArg", impact4.hgvsProtein());

        assertFalse(impact3.codingContext().IsFrameShift);
        assertFalse(impact3.hasEffect(FRAMESHIFT));
        assertTrue(impact3.phasedFrameshift());
        assertTrue(impact3.hasEffect(PHASED_INFRAME_INSERTION));
    }

    private VariantTransImpact classifyVariant(final VariantData var, final TranscriptData transData)
    {
        VariantTransImpact impact = mClassifier.classifyVariant(var, transData);
        var.addImpact(GENE_NAME_1, impact);
        return impact;
    }
}
