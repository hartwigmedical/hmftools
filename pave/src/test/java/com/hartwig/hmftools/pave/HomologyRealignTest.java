package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Codons.START_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.START_CODON;
import static com.hartwig.hmftools.common.codon.Codons.STOP_CODON_1;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INTRONIC;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SPLICE_ACCEPTOR;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createMockGenome;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateTestBases;
import static com.hartwig.hmftools.pave.ImpactTestUtils.getAminoAcidCodon;
import static com.hartwig.hmftools.pave.ImpactTestUtils.getAminoAcidsCodons;
import static com.hartwig.hmftools.pave.impact.PaveUtils.findVariantImpacts;
import static com.hartwig.hmftools.pave.impact.PaveUtils.createRightAlignedVariant;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import org.junit.Test;

public class HomologyRealignTest
{
    @Test
    public void testRealignment()
    {
        MockRefGenome refGenome = createMockGenome(201);
        String refBases = refGenome.RefGenomeMap.get(CHR_1);
        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        TranscriptData transDataPos = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] {20, 50, 80}, 10, 25, 65, false, "");

        GeneData geneData = createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, 20, 90);

        GeneDataCache geneDataCache = new GeneDataCache("", V37, null);
        geneDataCache.getEnsemblCache().getChrGeneDataMap().put(CHR_1, Lists.newArrayList(geneData));
        geneDataCache.getEnsemblCache().getTranscriptDataMap().put(GENE_ID_1, Lists.newArrayList(transDataPos));

        // intronic to splice
        int pos = 40; // intronic DEL
        String ref = refBases.substring(pos, pos + 2);
        String alt = refBases.substring(pos, pos + 1);

        VariantData var = new VariantData(CHR_1, pos, ref, alt);
        String altBases = ref.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 9);

        var.setRealignedVariant(createRightAlignedVariant(var, refGenome));
        assertTrue(var.realignedVariant() != null);

        findVariantImpacts(var, classifier, geneDataCache);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        assertEquals(1, var.getImpacts().get(GENE_NAME_1).size());

        VariantTransImpact impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertFalse(impact.realigned());
        assertEquals(INTRONIC, impact.topEffect());

        // now realigning from previous donor to intron
        pos = 33; // intronic DEL
        ref = refBases.substring(pos, pos + 2);
        alt = refBases.substring(pos, pos + 1);

        var = new VariantData(CHR_1, pos, ref, alt);
        altBases = ref.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 9);

        var.setRealignedVariant(createRightAlignedVariant(var, refGenome));
        assertTrue(var.realignedVariant() != null);

        findVariantImpacts(var, classifier, geneDataCache);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        assertEquals(1, var.getImpacts().get(GENE_NAME_1).size());

        impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertTrue(impact.realigned());
        assertEquals(INTRONIC, impact.topEffect());

        // now splice to exonic inframe INS
        pos = 48;
        ref = refBases.substring(pos, pos + 1);
        altBases = "AAA";
        alt = ref + altBases;

        var = new VariantData(CHR_1, pos, ref, alt);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 2);

        var.setRealignedVariant(createRightAlignedVariant(var, refGenome));
        assertTrue(var.realignedVariant() != null);

        findVariantImpacts(var, classifier, geneDataCache);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        assertEquals(1, var.getImpacts().get(GENE_NAME_1).size());

        impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertTrue(impact.realigned());
        assertEquals(INFRAME_INSERTION, impact.topEffect());
    }

    @Test
    public void testSpliceRealignment()
    {
        MockRefGenome refGenome = new MockRefGenome();

        int preGene = 10;
        int prePostCoding = 10;
        String refBases = generateTestBases(preGene);

        refBases += generateTestBases(prePostCoding);
        refBases += reverseComplementBases(STOP_CODON_1);
        refBases += "ATACCTGCT"; // Y-R-S going up, X (20-22), Y (23-25), R (26-28), S (29-31)
        refBases += "CTATAGAGCG"; // intron 32-41
        refBases += "CTTCTCCCT"; // K (42-44), E (45-47), R (48-50)
        refBases += reverseComplementBases(START_CODON); // M (51-52)
        refBases += generateTestBases(preGene);

        // post gene  3'UTR        X  Y  R  S intron       K  E  M 5'UTR
        // GATCGATCGA GATCGATCGA TTAATACCTGCT CTATAGAGCG CTTCTCCAT GATCGATCGA
        // 0          10         20           32

        TranscriptData negTransData =  createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, new int[] {10, 42}, 21,
                20, 52, false, "");

        refGenome.RefGenomeMap.put(CHR_1, refBases);

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // exonic to splice
        int pos = 26;
        String ref = refBases.substring(pos, pos + 6);
        String alt = refBases.substring(pos, pos + 1);

        VariantData var = new VariantData(CHR_1, pos, ref, alt);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, "CT", "C", 2);

        var.setRealignedVariant(createRightAlignedVariant(var, refGenome));
        assertTrue(var.realignedVariant() != null);

        GeneData geneData = createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, NEG_STRAND, negTransData.TransStart, negTransData.TransEnd);

        GeneDataCache geneDataCache = new GeneDataCache("", V37, null);
        geneDataCache.getEnsemblCache().getChrGeneDataMap().put(CHR_1, Lists.newArrayList(geneData));
        geneDataCache.getEnsemblCache().getTranscriptDataMap().put(GENE_ID_1, Lists.newArrayList(negTransData));

        findVariantImpacts(var, classifier, geneDataCache);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        assertEquals(1, var.getImpacts().get(GENE_NAME_1).size());

        VariantTransImpact impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertEquals(FRAMESHIFT, impact.topEffect());
        assertFalse(impact.realigned());

        VariantTransImpact raImpact = var.getRealignedImpact(GENE_NAME_1, impact);
        assertTrue(raImpact.hasEffect(SPLICE_ACCEPTOR));
    }

    @Test
    public void testSpliceVsFrameshiftRealignment()
    {
        MockRefGenome refGenome = new MockRefGenome();

        int prePostGene = 10;
        int prePostCoding = 10;
        String refBases = generateTestBases(prePostGene);

        refBases += generateTestBases(prePostCoding);

        String codingBases = getAminoAcidCodon(START_AMINO_ACID);
        codingBases += getAminoAcidsCodons("CFR", false);
        codingBases += getAminoAcidCodon('G', 2); // force to 'GGG'
        codingBases += getAminoAcidsCodons("D", false);
        codingBases += STOP_CODON_1;

        refBases += codingBases.substring(0, 10);
        refBases += generateTestBases(17) + "CAG"; // intron ending in splice acceptor motif with 'G'
        refBases += codingBases.substring(10);

        refBases += generateTestBases(prePostCoding);
        refBases += generateTestBases(prePostGene);

        // pre gene   5'UTR      M  C  F  R  intron               G  D  X  3'UTR
        // GATCGATCGA GATCGATCGA ATGTGTTTCA  GATCGATCGATCGATCGCAG GGGGGGACTAA GATCGATCGA GATCGATCGA
        // 0          10         20          30                   50

        TranscriptData posTransData =  createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] {10, 50}, 19,
                20, 60, false, "");

        refGenome.RefGenomeMap.put(CHR_1, refBases);

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // the initial variant should be splice acceptor, the second frameshift, so frameshift is selected by realignment priority rules
        int pos = 48;
        String ref = refBases.substring(pos, pos + 2);
        String alt = refBases.substring(pos, pos + 1);

        VariantData var = new VariantData(CHR_1, pos, ref, alt);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, "G", "G", 7);

        var.setRealignedVariant(createRightAlignedVariant(var, refGenome));
        assertTrue(var.realignedVariant() != null);

        GeneData geneData = createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, posTransData.TransStart, posTransData.TransEnd);

        GeneDataCache geneDataCache = new GeneDataCache("", V37, null);
        geneDataCache.getEnsemblCache().getChrGeneDataMap().put(CHR_1, Lists.newArrayList(geneData));
        geneDataCache.getEnsemblCache().getTranscriptDataMap().put(GENE_ID_1, Lists.newArrayList(posTransData));

        findVariantImpacts(var, classifier, geneDataCache);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        assertEquals(1, var.getImpacts().get(GENE_NAME_1).size());

        VariantTransImpact impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertEquals(FRAMESHIFT, impact.topEffect());
        assertTrue(impact.realigned());

        // now test an insertion of a base at A2 which keeps the splice acceptor motif the same but with a spare extra base before the exon
        pos = 48;
        ref = refBases.substring(pos, pos + 1);
        alt = "AG";

        var = new VariantData(CHR_1, pos, ref, alt);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, "G", "G", 7);

        var.setRealignedVariant(createRightAlignedVariant(var, refGenome));
        assertTrue(var.realignedVariant() != null);

        findVariantImpacts(var, classifier, geneDataCache);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        assertEquals(1, var.getImpacts().get(GENE_NAME_1).size());

        impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertEquals(FRAMESHIFT, impact.topEffect());
        assertTrue(impact.realigned());
    }
}
