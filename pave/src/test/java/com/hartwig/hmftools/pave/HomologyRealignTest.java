package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INTRONIC;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createMockGenome;
import static com.hartwig.hmftools.pave.PaveApplication.findVariantImpacts;
import static com.hartwig.hmftools.pave.PaveUtils.createRightAlignedVariant;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

import org.junit.Test;

public class HomologyRealignTest
{
    private MockRefGenome mRefGenome = createMockGenome(201);
    private String mRefBases = mRefGenome.RefGenomeMap.get(CHR_1);
    private ImpactClassifier mClassifier = new ImpactClassifier(mRefGenome);

    @Test
    public void testRealignment()
    {
        TranscriptData transDataPos = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] {20, 50, 80}, 10, 25, 65, false, "");

        GeneData geneData = createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, 20, 90);

        GeneDataCache geneDataCache = new GeneDataCache("", V37, null, false, false);
        geneDataCache.getEnsemblCache().getChrGeneDataMap().put(CHR_1, Lists.newArrayList(geneData));
        geneDataCache.getEnsemblCache().getTranscriptDataMap().put(GENE_ID_1, Lists.newArrayList(transDataPos));

        // intronic to splice
        int pos = 40; // intronic DEL
        String ref = mRefBases.substring(pos, pos + 2);
        String alt = mRefBases.substring(pos, pos + 1);

        VariantData var = new VariantData(CHR_1, pos, ref, alt);
        String altBases = ref.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 9);

        var.setRealignedVariant(createRightAlignedVariant(var, mRefGenome));
        assertTrue(var.realignedVariant() != null);

        findVariantImpacts(var, mClassifier, geneDataCache);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        assertEquals(1, var.getImpacts().get(GENE_NAME_1).size());

        VariantTransImpact impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertFalse(impact.realigned());
        assertEquals(INTRONIC, impact.topEffect());

        // now realigning from previous donor to intron
        pos = 33; // intronic DEL
        ref = mRefBases.substring(pos, pos + 2);
        alt = mRefBases.substring(pos, pos + 1);

        var = new VariantData(CHR_1, pos, ref, alt);
        altBases = ref.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 9);

        var.setRealignedVariant(createRightAlignedVariant(var, mRefGenome));
        assertTrue(var.realignedVariant() != null);

        findVariantImpacts(var, mClassifier, geneDataCache);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        assertEquals(1, var.getImpacts().get(GENE_NAME_1).size());

        impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertTrue(impact.realigned());
        assertEquals(INTRONIC, impact.topEffect());

        // now splice to exonic inframe INS
        pos = 48;
        ref = mRefBases.substring(pos, pos + 1);
        altBases = "AAA";
        alt = ref + altBases;

        var = new VariantData(CHR_1, pos, ref, alt);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 2);

        var.setRealignedVariant(createRightAlignedVariant(var, mRefGenome));
        assertTrue(var.realignedVariant() != null);

        findVariantImpacts(var, mClassifier, geneDataCache);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        assertEquals(1, var.getImpacts().get(GENE_NAME_1).size());

        impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertTrue(impact.realigned());
        assertEquals(INFRAME_INSERTION, impact.topEffect());
    }

}
