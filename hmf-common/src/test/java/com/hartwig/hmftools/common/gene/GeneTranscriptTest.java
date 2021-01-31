package com.hartwig.hmftools.common.gene;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.createBreakendTranscriptData;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.setAlternativeTranscriptPhasings;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.fusion.BreakendTransData.POST_CODING_PHASE;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.calcExonicCodingPhase;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.GeneTestUtils;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.common.fusion.CodingBaseData;

import org.junit.Assert;
import org.junit.Test;

public class GeneTranscriptTest
{
    private static final int TRANS_ID_1 = 1;
    private static final String TRANS_NAME_1 = "TRANS001";
    private static final String GENE_NAME_1 = "GENE_1";
    private static final String GENE_NAME_2 = "GENE_2";
    private static final String GENE_ID_1 = "GENE001";
    private static final String GENE_ID_2 = "GENE002";
    private static final String CHR_1 = "1";

    @Test
    public void testExonicCodingPhase()
    {
        ExonData exon = new ExonData(1, 10, 20, 1, PHASE_NONE, PHASE_NONE);

        assertEquals(PHASE_1, calcExonicCodingPhase(exon, 11, 19, POS_STRAND, 11));
        assertEquals(PHASE_1, calcExonicCodingPhase(exon, 11, 19, NEG_STRAND, 19));

        assertEquals(PHASE_1, calcExonicCodingPhase(exon, 11, 19, POS_STRAND, 14));
        assertEquals(PHASE_2, calcExonicCodingPhase(exon, 11, 19, POS_STRAND, 15));
        assertEquals(PHASE_2, calcExonicCodingPhase(exon, 11, 19, NEG_STRAND, 18));
        assertEquals(PHASE_0, calcExonicCodingPhase(exon, 11, 19, NEG_STRAND, 17));

        // special case of coding starting at the exon but with phase -1 specified
        assertEquals(PHASE_1, calcExonicCodingPhase(exon, 10, 25, POS_STRAND, 10));

        // invalid test since coding has started without phase being set for the exon data
        // assertEquals(PHASE_1, calcExonicCodingPhase(exon, 5, 20, NEG_STRAND, 20));

        // coding begins before the exon
        exon = new ExonData(1, 10, 20, 1, PHASE_2, PHASE_2);

        assertEquals(PHASE_0, calcExonicCodingPhase(exon, 5, 25, POS_STRAND, 10));
        assertEquals(PHASE_1, calcExonicCodingPhase(exon, 5, 25, POS_STRAND, 11));

        assertEquals(PHASE_0, calcExonicCodingPhase(exon, 5, 25, NEG_STRAND, 20));
        assertEquals(PHASE_1, calcExonicCodingPhase(exon, 5, 25, NEG_STRAND, 19));

        // coding begins at the exon boundary with phase specified
        assertEquals(PHASE_2, calcExonicCodingPhase(exon, 10, 25, POS_STRAND, 10));
        assertEquals(PHASE_0, calcExonicCodingPhase(exon, 10, 50, NEG_STRAND, 20));
    }

    @Test
    public void testCodingData()
    {
        Integer codingStart = null;
        Integer codingEnd = null;

        TranscriptData transData = new TranscriptData(
                TRANS_ID_1, TRANS_NAME_1, GENE_ID_1, false, POS_STRAND,
                10, 40, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        transData.exons().add(new ExonData(TRANS_ID_1, 10, 20, 1, -1, -1));
        transData.exons().add(new ExonData(TRANS_ID_1, 30, 40, 2, -1, -1));

        CodingBaseData cbData = calcCodingBases(transData, 25);

        assertEquals(0, cbData.TotalCodingBases);
        assertEquals(0, cbData.CodingBases);
        assertEquals(PHASE_NONE, cbData.Phase);

        codingStart = 35;
        codingEnd = 75;

        transData = new TranscriptData(
                TRANS_ID_1, TRANS_NAME_1, GENE_ID_1, false, POS_STRAND,
                10, 80, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        transData.exons().add(new ExonData(TRANS_ID_1, 10, 20, 1, PHASE_NONE, PHASE_NONE));
        transData.exons().add(new ExonData(TRANS_ID_1, 30, 40, 2, PHASE_NONE, PHASE_0));
        transData.exons().add(new ExonData(TRANS_ID_1, 50, 60, 2, PHASE_0, PHASE_2));
        transData.exons().add(new ExonData(TRANS_ID_1, 70, 80, 2, PHASE_2, PHASE_NONE));

        // pre-coding
        cbData = calcCodingBases(transData, 25);

        assertEquals(23, cbData.TotalCodingBases);
        assertEquals(0, cbData.CodingBases);
        assertEquals(PHASE_NONE, cbData.Phase);

        // post-coding
        cbData = calcCodingBases(transData, 76);

        assertEquals(cbData.TotalCodingBases, cbData.CodingBases);
        assertEquals(PHASE_NONE, cbData.Phase);

        // at first coding base
        cbData = calcCodingBases(transData, 35);

        assertEquals(1, cbData.CodingBases);
        assertEquals(PHASE_1, cbData.Phase);

        // at last coding base
        cbData = calcCodingBases(transData, 75);

        assertEquals(cbData.TotalCodingBases, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);

        // intronic
        cbData = calcCodingBases(transData, 45);

        assertEquals(6, cbData.CodingBases);
        assertEquals(PHASE_0, cbData.Phase);

        cbData = calcCodingBases(transData, 65);

        assertEquals(17, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);

        // test negative strand
        transData = new TranscriptData(
                TRANS_ID_1, TRANS_NAME_1, GENE_ID_1, false, NEG_STRAND,
                10, 80, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        transData.exons().add(new ExonData(TRANS_ID_1, 10, 20, 4, PHASE_NONE, PHASE_NONE));
        transData.exons().add(new ExonData(TRANS_ID_1, 30, 40, 3, PHASE_2, PHASE_NONE));
        transData.exons().add(new ExonData(TRANS_ID_1, 50, 60, 2, PHASE_0, PHASE_2));
        transData.exons().add(new ExonData(TRANS_ID_1, 70, 80, 1, PHASE_NONE, PHASE_0));

        cbData = calcCodingBases(transData, 25);

        assertEquals(23, cbData.TotalCodingBases);
        assertEquals(cbData.TotalCodingBases, cbData.CodingBases);
        assertEquals(PHASE_NONE, cbData.Phase);

        cbData = calcCodingBases(transData, 76);

        assertEquals(0, cbData.CodingBases);
        assertEquals(PHASE_NONE, cbData.Phase);

        cbData = calcCodingBases(transData, 35);

        assertEquals(cbData.TotalCodingBases, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);

        cbData = calcCodingBases(transData, 75);

        assertEquals(1, cbData.CodingBases);
        assertEquals(PHASE_1, cbData.Phase);

        cbData = calcCodingBases(transData, 45);

        assertEquals(17, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);

        cbData = calcCodingBases(transData, 65);

        assertEquals(6, cbData.CodingBases);
        assertEquals(PHASE_0, cbData.Phase);

        // test coding start on the first base of the first exon with specific phasing
        codingStart = 10;
        codingEnd = 50;
        transData = new TranscriptData(
                TRANS_ID_1, TRANS_NAME_1, GENE_ID_1, false, POS_STRAND,
                10, 80, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        transData.exons().add(new ExonData(TRANS_ID_1, 10, 20, 1, PHASE_2, PHASE_0));
        transData.exons().add(new ExonData(TRANS_ID_1, 30, 40, 2, PHASE_0, PHASE_1));
        transData.exons().add(new ExonData(TRANS_ID_1, 50, 60, 3, PHASE_1, PHASE_NONE));

        cbData = calcCodingBases(transData, 10);

        assertEquals(1, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);

        cbData = calcCodingBases(transData, 20);

        assertEquals(11, cbData.CodingBases);
        assertEquals(PHASE_0, cbData.Phase);

        // and on the negative strand
        codingStart = 30;
        codingEnd = 60;
        transData = new TranscriptData(
                TRANS_ID_1, TRANS_NAME_1, GENE_ID_1, false, NEG_STRAND,
                10, 80, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        transData.exons().add(new ExonData(TRANS_ID_1, 10, 20, 3, PHASE_NONE, PHASE_NONE));
        transData.exons().add(new ExonData(TRANS_ID_1, 30, 40, 2, PHASE_0, PHASE_1));
        transData.exons().add(new ExonData(TRANS_ID_1, 50, 60, 1, PHASE_2, PHASE_0));

        cbData = calcCodingBases(transData, 60);

        assertEquals(1, cbData.CodingBases);
        assertEquals(PHASE_2, cbData.Phase);

        cbData = calcCodingBases(transData, 50);

        assertEquals(11, cbData.CodingBases);
        assertEquals(PHASE_0, cbData.Phase);

        cbData = calcCodingBases(transData, 40);

        assertEquals(12, cbData.CodingBases);
        assertEquals(PHASE_1, cbData.Phase);
    }

    @Test
    public void testTranscriptDataCreation()
    {
        int[] exonStarts = {10, 30, 50};
        Integer codingStart = 15;
        Integer codingEnd = 55;

        TranscriptData transData = GeneTestUtils.createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, BIOTYPE_PROTEIN_CODING);

        assertEquals(3, transData.exons().size());
        assertEquals(10, transData.exons().get(0).Start);
        assertEquals(20, transData.exons().get(0).End);
        assertEquals(3, transData.exons().get(2).Rank);

        assertEquals(PHASE_NONE, transData.exons().get(0).PhaseStart);
        assertEquals(PHASE_0, transData.exons().get(0).PhaseEnd);
        assertEquals(PHASE_0, transData.exons().get(1).PhaseStart);
        assertEquals(PHASE_2, transData.exons().get(1).PhaseEnd);
        assertEquals(PHASE_2, transData.exons().get(2).PhaseStart);
        assertEquals(PHASE_NONE, transData.exons().get(2).PhaseEnd);

        transData = GeneTestUtils.createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, BIOTYPE_PROTEIN_CODING);

        assertEquals(3, transData.exons().size());
        assertEquals(50, transData.exons().get(2).Start);
        assertEquals(60, transData.exons().get(2).End);
        assertEquals(3, transData.exons().get(0).Rank);

        assertEquals(PHASE_NONE, transData.exons().get(2).PhaseStart);
        assertEquals(PHASE_0, transData.exons().get(2).PhaseEnd);
        assertEquals(PHASE_0, transData.exons().get(1).PhaseStart);
        assertEquals(PHASE_2, transData.exons().get(1).PhaseEnd);
        assertEquals(PHASE_2, transData.exons().get(0).PhaseStart);
        assertEquals(PHASE_NONE, transData.exons().get(0).PhaseEnd);

        // coding start on the first base of an exon, is given the phase before the base
        transData = GeneTestUtils.createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, 10, codingEnd, false, BIOTYPE_PROTEIN_CODING);

        assertEquals(PHASE_0, transData.exons().get(0).PhaseStart);
        assertEquals(PHASE_2, transData.exons().get(0).PhaseEnd);

        transData = GeneTestUtils.createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, 15, 60, false, BIOTYPE_PROTEIN_CODING);

        assertEquals(PHASE_0, transData.exons().get(2).PhaseStart);
        assertEquals(PHASE_2, transData.exons().get(2).PhaseEnd);
    }

    @Test
    public void testAlternativeBreakendTranscriptPhasing()
    {
        EnsemblDataCache geneTransCache = createGeneDataCache();

        // first a gene on the forward strand
        String geneName = "GENE1";
        String geneId = "ENSG0001";

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(GeneTestUtils.createEnsemblGeneData(geneId, geneName, CHR_1, POS_STRAND, 100, 1000));
        GeneTestUtils.addGeneData(geneTransCache, CHR_1, geneList);

        BreakendGeneData genePosStrand = GeneTestUtils.createGeneAnnotation(0, true, geneName, geneId, POS_STRAND, CHR_1, 0, POS_ORIENT);

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

        BreakendGeneData geneNegStrand = GeneTestUtils.createGeneAnnotation(0, true, geneName, geneId, NEG_STRAND, CHR_1, 0, POS_ORIENT);

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
        exonsSkipped = trans.getAlternativePhasing().get(PHASE_NONE);
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
        BreakendGeneData genePosStrand = GeneTestUtils.createGeneAnnotation(0, true, GENE_NAME_1, GENE_ID_1, POS_STRAND, CHR_1, 0, POS_ORIENT);

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
        BreakendGeneData geneNegStrand = GeneTestUtils.createGeneAnnotation(0, true, GENE_NAME_2, GENE_ID_2, NEG_STRAND, CHR_1, 0, POS_ORIENT);

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
