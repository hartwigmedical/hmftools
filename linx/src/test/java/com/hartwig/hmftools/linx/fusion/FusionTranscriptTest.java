package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.gene.BreakendGenePrep.createBreakendTranscriptData;
import static com.hartwig.hmftools.linx.gene.BreakendGenePrep.setAlternativeTranscriptPhasings;
import static com.hartwig.hmftools.linx.gene.BreakendTransData.POST_CODING_PHASE;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.createGeneAnnotation;

import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.test.GeneTestUtils;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.linx.gene.BreakendTransData;

import org.junit.Assert;
import org.junit.Test;

public class FusionTranscriptTest
{
    private static final String GENE_NAME_1 = "GENE_1";
    private static final String GENE_NAME_2 = "GENE_2";
    private static final String GENE_ID_1 = "GENE001";
    private static final String GENE_ID_2 = "GENE002";
    private static final String CHR_1 = "1";

    @Test
    public void testAlternativeBreakendTranscriptPhasing()
    {
        EnsemblDataCache geneTransCache = createGeneDataCache();

        // first a gene on the forward strand
        String geneName = "GENE1";
        String geneId = "ENSG0001";

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(GeneTestUtils.createEnsemblGeneData(geneId, geneName, CHR_1, POS_STRAND, 100, 1000));
        GeneTestUtils.addGeneData(geneTransCache, CHR_1, geneList);

        BreakendGeneData genePosStrand = createGeneAnnotation(0, true, geneName, geneId, POS_STRAND, CHR_1, 0, POS_ORIENT);

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;
        int[] exonStarts = new int[]{100, 300, 500, 700, 900};

        Integer codingStart = 349;
        Integer codingEnd = 950;
        TranscriptData transData = createTransExons(geneId, transId++, POS_STRAND, exonStarts, 100, codingStart, codingEnd, false, "");
        transDataList.add(transData);

        GeneTestUtils.addTransExonData(geneTransCache, geneId, transDataList);

        int position = 250;
        genePosStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        BreakendTransData trans = createBreakendTranscriptData(transData, position, genePosStrand);

        Assert.assertEquals(5, trans.exonCount());
        Assert.assertEquals(1, trans.ExonUpstream);
        Assert.assertEquals(2, trans.ExonDownstream);
        Assert.assertEquals(PHASE_NONE, trans.Phase);

        // test caching of upstream phasings for exon-skipping fusion logic
        setAlternativeTranscriptPhasings(trans, transData.exons(), position, POS_ORIENT);
        Assert.assertEquals(0, trans.getAlternativePhasing().size());

        // and test as a downstream gene
        setAlternativeTranscriptPhasings(trans, transData.exons(), position, NEG_ORIENT);
        Assert.assertEquals(3, trans.getAlternativePhasing().size());
        Integer exonsSkipped = trans.getAlternativePhasing().get(PHASE_1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_2);
        assertTrue(exonsSkipped != null && exonsSkipped == 3);
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_0);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);

        position = 450;
        genePosStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        trans = createBreakendTranscriptData(transData, position, genePosStrand);

        Assert.assertEquals(2, trans.ExonUpstream);
        Assert.assertEquals(3, trans.ExonDownstream);
        Assert.assertEquals(PHASE_1, trans.Phase);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, POS_ORIENT);
        Assert.assertEquals(1, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(-1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, NEG_ORIENT);
        Assert.assertEquals(2, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_2);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_0);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);

        position = 650;
        genePosStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        trans = createBreakendTranscriptData(transData, position, genePosStrand);

        Assert.assertEquals(3, trans.ExonUpstream);
        Assert.assertEquals(4, trans.ExonDownstream);
        Assert.assertEquals(PHASE_0, trans.Phase);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, POS_ORIENT);
        Assert.assertEquals(2, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(-1);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, NEG_ORIENT);
        Assert.assertEquals(1, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_2);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);

        // then a gene on the reverse strand
        geneName = "GENE2";
        geneId = "ENSG0002";

        geneList.add(GeneTestUtils.createEnsemblGeneData(geneId, geneName, CHR_1, POS_STRAND, 100, 1000));
        GeneTestUtils.addGeneData(geneTransCache, CHR_1, geneList);

        BreakendGeneData geneNegStrand = createGeneAnnotation(0, true, geneName, geneId, NEG_STRAND, CHR_1, 0, POS_ORIENT);

        transDataList = Lists.newArrayList();

        transId = 2;

        exonStarts = new int[]{100, 300, 500, 700, 900};
        transData = createTransExons(geneId, transId++, NEG_STRAND, exonStarts, 100, codingStart, codingEnd, false, "");

        transDataList.add(transData);

        GeneTestUtils.addTransExonData(geneTransCache, geneId, transDataList);

        position = 850;
        geneNegStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        trans = createBreakendTranscriptData(transData, position, geneNegStrand);

        Assert.assertEquals(5, trans.exonCount());
        Assert.assertEquals(1, trans.ExonUpstream);
        Assert.assertEquals(2, trans.ExonDownstream);
        Assert.assertEquals(PHASE_0, trans.Phase);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, NEG_ORIENT);
        Assert.assertEquals(0, trans.getAlternativePhasing().size());

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, POS_ORIENT);
        Assert.assertEquals(3, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_1);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);
        exonsSkipped = trans.getAlternativePhasing().get(POST_CODING_PHASE);
        assertTrue(exonsSkipped != null && exonsSkipped == 3);
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_2);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);

        position = 250;
        geneNegStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        trans = createBreakendTranscriptData(transData, position, geneNegStrand);

        Assert.assertEquals(4, trans.ExonUpstream);
        Assert.assertEquals(5, trans.ExonDownstream);
        Assert.assertEquals(POST_CODING_PHASE, trans.Phase);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, NEG_ORIENT);
        Assert.assertEquals(3, trans.getAlternativePhasing().size());
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_1);
        assertTrue(exonsSkipped != null && exonsSkipped == 1);
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_0);
        assertTrue(exonsSkipped != null && exonsSkipped == 3);
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_2);
        assertTrue(exonsSkipped != null && exonsSkipped == 2);

        setAlternativeTranscriptPhasings(trans, transData.exons(), position, POS_ORIENT);
        Assert.assertEquals(0, trans.getAlternativePhasing().size());
    }

    @Test
    public void testBreakendTranscriptCoding()
    {
        // first a gene on the forward strand
        BreakendGeneData genePosStrand = createGeneAnnotation(0, true, GENE_NAME_1, GENE_ID_1, POS_STRAND, CHR_1, 0, POS_ORIENT);

        int transId = 1;

        int[] exonStarts = new int[]{100, 200, 300, 400, 500};

        // coding taking up exactly the first exon
        Integer codingStart = 100;
        Integer codingEnd = 110;

        TranscriptData transData = createTransExons(
                GENE_ID_1, transId++, POS_STRAND, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        int position = 150;
        genePosStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        BreakendTransData trans = createBreakendTranscriptData(transData, position, genePosStrand);

        Assert.assertEquals(5, trans.exonCount());
        Assert.assertEquals(1, trans.ExonUpstream);
        Assert.assertEquals(2, trans.ExonDownstream);
        Assert.assertEquals(UTR_3P, trans.codingType());
        Assert.assertEquals(-2, trans.Phase);
        Assert.assertEquals(-2, trans.Phase);
        Assert.assertEquals(8, trans.CodingBases); // stop codon is taken out
        Assert.assertEquals(8, trans.TotalCodingBases);

        //
        codingStart = 105;
        codingEnd = 405;

        transData = createTransExons(
                GENE_ID_1, transId++, POS_STRAND, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        position = 350;
        genePosStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        trans = createBreakendTranscriptData(transData, position, genePosStrand);

        Assert.assertEquals(3, trans.ExonUpstream);
        Assert.assertEquals(4, trans.ExonDownstream);
        Assert.assertEquals(CODING, trans.codingType());
        Assert.assertEquals(PHASE_1, trans.Phase);
        Assert.assertEquals(28, trans.CodingBases);
        Assert.assertEquals(31, trans.TotalCodingBases);

        // test the reverse strand
        BreakendGeneData geneNegStrand = createGeneAnnotation(0, true, GENE_NAME_2, GENE_ID_2, NEG_STRAND, CHR_1, 0, POS_ORIENT);

        // coding taking up exactly the first exon
        codingStart = 500;
        codingEnd = 510;

        transData = createTransExons(
                GENE_ID_2, transId++, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        position = 450;
        geneNegStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        trans = createBreakendTranscriptData(transData, position, geneNegStrand);

        Assert.assertEquals(5, trans.exonCount());
        Assert.assertEquals(1, trans.ExonUpstream);
        Assert.assertEquals(2, trans.ExonDownstream);
        Assert.assertEquals(UTR_3P, trans.codingType());
        Assert.assertEquals(-2, trans.Phase);
        Assert.assertEquals(8, trans.CodingBases);
        Assert.assertEquals(8, trans.TotalCodingBases);

        codingStart = 205;
        codingEnd = 505;

        transData = createTransExons(
                GENE_ID_2, transId++, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        position = 250;
        geneNegStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        trans = createBreakendTranscriptData(transData, position, geneNegStrand);

        Assert.assertEquals(3, trans.ExonUpstream);
        Assert.assertEquals(4, trans.ExonDownstream);
        Assert.assertEquals(CODING, trans.codingType());
        Assert.assertEquals(PHASE_1, trans.Phase);
        Assert.assertEquals(28, trans.CodingBases);
        Assert.assertEquals(31, trans.TotalCodingBases);

        // test coding starting on the first base of the second exon for a downstream transcript
        codingStart = 200;
        codingEnd = 405;

        transData = createTransExons(
                GENE_ID_1, transId++, POS_STRAND, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        position = 50;
        genePosStrand.setPositionalData(CHR_1, position, NEG_ORIENT);
        trans = createBreakendTranscriptData(transData, position, genePosStrand);

        Assert.assertEquals(0, trans.ExonUpstream);
        Assert.assertEquals(2, trans.ExonDownstream);
        Assert.assertEquals(UTR_5P, trans.codingType());
        Assert.assertEquals(PHASE_NONE, trans.Phase);

        // coding starting within the first exon
        codingStart = 105;
        codingEnd = 405;

        transData = createTransExons(
                GENE_ID_1, transId++, POS_STRAND, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        trans = createBreakendTranscriptData(transData, position, genePosStrand);

        Assert.assertEquals(0, trans.ExonUpstream);
        Assert.assertEquals(2, trans.ExonDownstream);
        Assert.assertEquals(CODING, trans.codingType());
        Assert.assertEquals(PHASE_0, trans.Phase);

        // test coding starting on a subsequent exon for a downstream transcript
        codingStart = 105;
        codingEnd = 310;

        transData = createTransExons(
                GENE_ID_2, transId++, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        position = 350;
        geneNegStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        trans = createBreakendTranscriptData(transData, position, geneNegStrand);

        Assert.assertEquals(2, trans.ExonUpstream);
        Assert.assertEquals(3, trans.ExonDownstream);
        Assert.assertEquals(UTR_5P, trans.codingType());
        Assert.assertEquals(PHASE_NONE, trans.Phase);

        codingStart = 105;
        codingEnd = 405;

        transData = createTransExons(
                GENE_ID_2, transId++, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, true, BIOTYPE_PROTEIN_CODING);

        position = 350;
        geneNegStrand.setPositionalData(CHR_1, position, POS_ORIENT);
        trans = createBreakendTranscriptData(transData, position, geneNegStrand);

        Assert.assertEquals(2, trans.ExonUpstream);
        Assert.assertEquals(3, trans.ExonDownstream);
        Assert.assertEquals(CODING, trans.codingType());
        Assert.assertEquals(PHASE_0, trans.Phase);
    }
}
