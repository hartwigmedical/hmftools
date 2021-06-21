package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneAnnotation;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.getCodingBases;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.BreakendTransData.POST_CODING_PHASE;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.REQUIRED_BIOTYPES;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.checkFusionLogic;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_ID_2;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_NAME_2;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.TRANS_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.TRANS_2;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.createTranscript;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

import junit.framework.TestCase;

public class FusionRulesTest
{
    @Test
    public void testFusionTypes()
    {
        String geneName = "GENE1";
        String geneId = "ENSG0001";

        // SV breakend positions won't impact fusion determination since transcripts are created manually
        BreakendGeneData gene1 = createGeneAnnotation(0, true, geneName, geneId, POS_STRAND, CHR_1, 150, POS_ORIENT);

        // one on the negative strand
        String geneName2 = "GENE2";
        String geneId2 = "ENSG0003";

        BreakendGeneData gene2 = createGeneAnnotation(0, false, geneName2, geneId2, NEG_STRAND, CHR_1, 150, POS_ORIENT);

        int transId1 = 1;

        /* invalid fusions:
            - up or down post-coding
            - down non-coding
            - up promoter
            - down single exon
            - up not disruptive
            - up 5'UTR to down coding
            - up coding to down not-coding
            - up coding exonic to not down coding exonic
            - up non-coding to down coding
            - unmatched phasing - exact or not
         */

        // non-coding combos
        Integer codingStart = null;
        Integer codingEnd = null;

        BreakendTransData trans1 = createTranscript(
                gene1, transId1, false, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_NONE, 0, getCodingBases(codingStart, codingEnd));

        int transId2 = 2;

        codingStart = 100;
        codingEnd = 200;

        BreakendTransData trans2 = createTranscript(
                gene2, transId2, false, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_NONE, 10, getCodingBases(codingStart, codingEnd));

        FusionParameters params = new FusionParameters();
        params.AllowExonSkipping = false;
        params.RequirePhaseMatch = true;

        // up non-coding
        assertTrue(trans1.nonCoding());
        assertTrue(trans2.isCoding());
        assertNull(checkFusionLogic(trans1, trans2, params));

        codingStart = null;
        codingEnd = null;

        trans2 = createTranscript(
                gene2, transId2, false, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_NONE, 10, getCodingBases(codingStart, codingEnd));

        codingStart = 100;
        codingEnd = 200;
        trans1 = createTranscript(
                gene2, transId2, false, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_NONE, 10, getCodingBases(codingStart, codingEnd));

        // down non-coding
        assertTrue(trans1.isCoding());
        assertTrue(trans2.nonCoding());
        assertNull(checkFusionLogic(trans1, trans2, params));

        // up / down post coding - upstream can be post-coding if can skip exons to form a phase-matched fusion
        trans2 = createTranscript(
                gene2, transId2, false, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_0, 10, getCodingBases(codingStart, codingEnd));

        gene1.setPositionalData(CHR_1, 260, POS_ORIENT);
        trans1 = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_NONE, getCodingBases(codingStart, codingEnd), getCodingBases(codingStart, codingEnd));

        assertTrue(trans1.postCoding());
        assertTrue(trans2.isCoding());
        assertNull(checkFusionLogic(trans1, trans2, params));

        // up promotor
        trans1 = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                0, 1, PHASE_NONE, 0, getCodingBases(codingStart, codingEnd));

        assertTrue(trans1.isPromoter());
        assertNull(checkFusionLogic(trans1, trans2, params));

        // down single exon
        gene1.setPositionalData(CHR_1, 40, POS_ORIENT);
        trans1 = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_NONE, 0, getCodingBases(codingStart, codingEnd));

        trans2 = createTranscript(
                gene2, transId2, true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_NONE, 10, getCodingBases(codingStart, codingEnd));

        assertTrue(trans1.preCoding());
        assertNull(checkFusionLogic(trans1, trans2, params));

        // up not disruptive
        trans1.setIsDisruptive(false);

        trans2 = createTranscript(
                gene2, transId2,  true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_NONE, 10, getCodingBases(codingStart, codingEnd));

        assertNull(checkFusionLogic(trans1, trans2, params));

        // up 5'UTR to down coding
        trans1.setIsDisruptive(true);

        // up coding exonic to not down coding exonic
        trans1 = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                2, 2, PHASE_NONE, 0, getCodingBases(codingStart, codingEnd));

        trans2 = createTranscript(
                gene2, transId2, false, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_NONE, 0, getCodingBases(codingStart, codingEnd));

        assertTrue(trans1.isExonic());
        assertTrue(trans2.isIntronic());
        assertNull(checkFusionLogic(trans1, trans2, params));

        // unmatched phasing intronic
        gene1.setPositionalData(CHR_1, 110, POS_ORIENT);
        trans1 = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_NONE, 10, getCodingBases(codingStart, codingEnd));

        trans2 = createTranscript(
                gene2, transId2,  true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_0, 10, getCodingBases(codingStart, codingEnd));

        assertTrue(trans1.isCoding());
        assertTrue(trans2.isCoding());
        assertTrue(trans1.isIntronic());
        assertTrue(trans2.isIntronic());
        assertNull(checkFusionLogic(trans1, trans2, params));

        // pre-coding to pre-coding same gene
        gene1.setPositionalData(CHR_1, 40, POS_ORIENT);
        trans1 = createTranscript(
                gene1, transId1, true, 50, 450, codingStart, codingEnd, "",
                1, 2, PHASE_NONE, 0, getCodingBases(codingStart, codingEnd));

        trans2 = createTranscript(
                gene1, transId1, true, 50, 550, codingStart, codingEnd, "",
                4, 5, PHASE_NONE, 0, getCodingBases(codingStart, codingEnd));

        assertTrue(trans1.preCoding());
        assertTrue(trans2.preCoding());
        assertNull(checkFusionLogic(trans1, trans2, params));

        // exon to exon phasing - to be phased, the upstream phase needs to be 1 less than the downstream phase

        // upstream coding bases - 10, phase = 0, will need the downstream phase to be 1

        // unmatched phasing exon to exon
        gene1.setPositionalData(CHR_1, 110, POS_ORIENT);
        trans1 = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                2, 2, PHASE_1, PHASE_2, 10, getCodingBases(codingStart, codingEnd));

        gene2.setPositionalData(CHR_1, 190, POS_ORIENT);
        trans2 = createTranscript(
                gene2, transId2, true, 50, 250, codingStart, codingEnd, "",
                2, 2, PHASE_0, PHASE_2, 10, getCodingBases(codingStart, codingEnd));

        assertTrue(trans1.isExonic());
        assertTrue(trans2.isExonic());
        assertNull(checkFusionLogic(trans1, trans2, params));

        // also invalid exonic fusion
        trans1 = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                2, 2, PHASE_2, 11, getCodingBases(codingStart, codingEnd));

        assertNull(checkFusionLogic(trans1, trans2, params));

        // valid exonic fusion
        trans2 = createTranscript(
                gene2, transId2, true, 50, 250, codingStart, codingEnd, "",
                2, 2, PHASE_0, 12, getCodingBases(codingStart, codingEnd));

        assertNotNull(checkFusionLogic(trans1, trans2, params));

        // exon to exon but both pre-coding
        gene1.setPositionalData(CHR_1, 90, POS_ORIENT);
        trans1 = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                2, 2, PHASE_NONE, 0, getCodingBases(codingStart, codingEnd));

        gene2.setPositionalData(CHR_1, 210, POS_ORIENT);
        trans2 = createTranscript(
                gene2, transId2, true, 50, 250, codingStart, codingEnd, "",
                2, 2, PHASE_NONE, 0, getCodingBases(codingStart, codingEnd));

        assertTrue(trans1.isExonic());
        assertTrue(trans2.isExonic());
        assertTrue(trans1.preCoding());
        assertTrue(trans2.preCoding());
        assertNotNull(checkFusionLogic(trans1, trans2, params));

        // valid intronic fusion
        trans1 = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_1, 10, getCodingBases(codingStart, codingEnd));

        trans2 = createTranscript(
                gene2, transId2, true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_1, 10, getCodingBases(codingStart, codingEnd));

        assertNotNull(checkFusionLogic(trans1, trans2, params));

        FusionFinder fusionFinder = new FusionFinder(null, null);
        fusionFinder.getKnownFusionCache().addData(new KnownFusionData(KNOWN_PAIR, gene1.GeneName, gene2.GeneName, "", ""));
        gene1.addTranscript(trans1);
        gene2.addTranscript(trans2);
        final List<GeneFusion> fusions = fusionFinder.findFusions(Lists.newArrayList(gene1), Lists.newArrayList(gene2), params);
        assertEquals(1, fusions.size());
        assertTrue(fusions.get(0).phaseMatched());
    }

    @Test
    public void testAlternatePhasings()
    {
        String geneName = "GENE1";
        String geneId = "ENSG0001";

        // SV breakend positions won't impact fusion determination since transcripts are created manually
        BreakendGeneData gene1 = createGeneAnnotation(0, true, geneName, geneId, POS_STRAND, CHR_1, 150, POS_ORIENT);

        // one on the negative strand
        String geneName2 = "GENE2";
        String geneId2 = "ENSG0003";

        BreakendGeneData gene2 = createGeneAnnotation(0, false, geneName2, geneId2, NEG_STRAND, CHR_1, 150, POS_ORIENT);

        int transId1 = 1;

        // non-coding combos
        Integer codingStart = 100;
        Integer codingEnd = 200;

        BreakendTransData transUp = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_1, 10, getCodingBases(codingStart, codingEnd));

        int transId2 = 2;

        BreakendTransData transDown = createTranscript(
                gene2, transId2, false, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_0, 10, getCodingBases(codingStart, codingEnd));

        FusionParameters params = new FusionParameters();
        params.AllowExonSkipping = true;
        params.RequirePhaseMatch = false;

        // up non-coding
        assertTrue(transUp.isCoding());
        assertTrue(transDown.isCoding());
        GeneFusion fusion = checkFusionLogic(transUp, transDown, params);

        assertNotNull(fusion);
        assertFalse(fusion.phaseMatched());

        Map<Integer,Integer> altPhasingsUp = transUp.getAlternativePhasing();
        altPhasingsUp.put(0, 1);

        fusion = checkFusionLogic(transUp, transDown, params);

        assertNotNull(fusion);
        assertTrue(fusion.phaseMatched());
        assertEquals(fusion.getExonsSkipped(true), 1);
        assertEquals(fusion.getExonsSkipped(false), 0);

        altPhasingsUp.clear();
        Map<Integer,Integer> altPhasingsDown = transDown.getAlternativePhasing();
        altPhasingsDown.put(1, 1);

        fusion = checkFusionLogic(transUp, transDown, params);

        assertNotNull(fusion);
        assertTrue(fusion.phaseMatched());
        assertEquals(0, fusion.getExonsSkipped(true));
        assertEquals(1, fusion.getExonsSkipped(false));

        // check 5' gene fusing from the 3'UTR region
        gene1.setPositionalData(CHR_1, 210, POS_ORIENT);
        transUp = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                6, 7, PHASE_NONE, 100, 100);

        assertTrue(transUp.postCoding());
        assertEquals(transUp.Phase, POST_CODING_PHASE);

        altPhasingsUp = transUp.getAlternativePhasing();

        altPhasingsUp.clear();
        altPhasingsUp.put(0, 3);

        fusion = checkFusionLogic(transUp, transDown, params);

        assertNotNull(fusion);
        assertTrue( fusion.phaseMatched());
        assertEquals(fusion.getExonsSkipped(true), 3);
        assertEquals(fusion.getExonsSkipped(false), 0);

        // check a fusion requiring skipping on both transcripts
        transUp = createTranscript(
                gene1, transId1, true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_1, 10, getCodingBases(codingStart, codingEnd));;

        transDown = createTranscript(
                gene2, transId2, true, 50, 250, codingStart, codingEnd, "",
                2, 3, PHASE_0, 10, getCodingBases(codingStart, codingEnd));;

        params.AllowExonSkipping = true;

        fusion = checkFusionLogic(transUp, transDown, params);

        assertNotNull(fusion);
        assertFalse(fusion.phaseMatched());

        altPhasingsUp = transUp.getAlternativePhasing();
        altPhasingsDown = transDown.getAlternativePhasing();
        altPhasingsUp.put(2, 2);
        altPhasingsDown.put(2, 3);

        fusion = checkFusionLogic(transUp, transDown, params);

        assertNotNull(fusion);
        assertTrue(fusion.phaseMatched());
        assertEquals(2, fusion.getExonsSkipped(true));
        assertEquals(3, fusion.getExonsSkipped(false));
    }

    @Test
    public void testExonicFusions()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();
        tester.initialiseFusions(geneTransCache);

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, 100, 1000));

        int[] exonStarts = new int[]{100, 300, 500, 700, 900};

        TranscriptData transData = createTransExons(GENE_ID_1, TRANS_1, POS_STRAND, exonStarts,  100,
                150, 950, true, BIOTYPE_PROTEIN_CODING);

        addTransExonData(geneTransCache, GENE_ID_1, Lists.newArrayList(transData));

        geneList.add(createEnsemblGeneData(GENE_ID_2, GENE_NAME_2, CHR_1, POS_STRAND, 10100, 11000));

        addGeneData(geneTransCache, CHR_1, geneList);

        exonStarts = new int[]{10100, 10300, 10500, 10700, 10900};

        transData = createTransExons(GENE_ID_2, TRANS_2, POS_STRAND, exonStarts, 100,
                10150, 10950, true, BIOTYPE_PROTEIN_CODING);

        addTransExonData(geneTransCache, GENE_ID_2, Lists.newArrayList(transData));

        // upstream exonic to downstream intronic - will use the preceding exon and exon's end phase to fuse
        List<BreakendGeneData> upGenes = Lists.newArrayList();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, CHR_1, 701, POS_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(CHR_1, 701, POS_ORIENT);

        List<BreakendGeneData> downGenes = Lists.newArrayList();
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, CHR_1, 10650, NEG_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(CHR_1, 10650, NEG_ORIENT);

        FusionParameters params = new FusionParameters();
        params.RequirePhaseMatch = false;
        params.AllowExonSkipping = true;

        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        assertEquals(1, fusions.size());
        GeneFusion fusion = fusions.get(0);

        // the selected fusion is the longest for coding bases and without any exon skipping
        assertTrue(fusion.phaseMatched());
        assertEquals(0, fusion.transcripts()[FS_UP].ExonicBasePhase);
        assertEquals(1, fusion.transcripts()[FS_UP].Phase);
        assertEquals(EXONIC, fusion.transcripts()[FS_UP].regionType());
        assertEquals(4, fusion.transcripts()[FS_UP].ExonUpstream);
        assertEquals(3, fusion.transcripts()[FS_DOWN].ExonUpstream);
        assertEquals(3, fusion.getFusedExon(true));
        assertEquals(0, fusion.getExonsSkipped(true));
        assertEquals(0, fusion.getExonsSkipped(false));

        // now an exonic fusion
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, CHR_1, 701, POS_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(CHR_1, 701, POS_ORIENT);

        downGenes.clear();
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, CHR_1, 10705, NEG_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(CHR_1, 10705, NEG_ORIENT);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        assertEquals(1, fusions.size());
        fusion = fusions.get(0);

        // the selected fusion is the longest for coding bases and without any exon skipping
        assertTrue(fusion.phaseMatched());
        assertEquals(0, fusion.transcripts()[FS_UP].ExonicBasePhase);
        assertEquals(1, fusion.transcripts()[FS_DOWN].ExonicBasePhase);
        assertEquals(EXONIC, fusion.transcripts()[FS_UP].regionType());
        assertEquals(EXONIC, fusion.transcripts()[FS_DOWN].regionType());
        assertEquals(4, fusion.transcripts()[FS_UP].ExonUpstream);
        assertEquals(4, fusion.transcripts()[FS_DOWN].ExonUpstream);
        assertEquals(4, fusion.getFusedExon(true));
        assertEquals(4, fusion.getFusedExon(false));
    }
}
