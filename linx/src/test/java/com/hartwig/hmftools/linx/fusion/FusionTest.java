package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneAnnotation;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_PROCESSED_TRANS;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.analysis.VariantPrep.setSvGeneData;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.findTopPriorityFusion;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.createTranscript;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.linx.SampleAnalyser;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class FusionTest
{
    @Test
    public void testReportableFusionComparison()
    {
        // test the selection and prioritisation logic for reportable fusions
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        PRE_GENE_PROMOTOR_DISTANCE = 100;

        // first a gene on the forward strand
        String geneName = "GENE1";
        String geneId1 = "ENSG0001";

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, geneName, CHR_1, POS_STRAND, 100, 1000));

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;

        int[] exonStarts = new int[]{100, 300, 500, 700, 900};

        int codingStart = 300;
        int codingEnd = 998;
        TranscriptData transData = createTransExons(geneId1, transId++, POS_STRAND, exonStarts,  99, codingStart, codingEnd, true, "");
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId1, transDataList);

        geneName = "GENE2";
        String geneId2 = "ENSG0002";

        geneList.add(createEnsemblGeneData(geneId2, geneName, CHR_1, 1, 10000, 12000));

        addGeneData(geneTransCache, CHR_1, geneList);

        transDataList = Lists.newArrayList();

        exonStarts = new int[]{10100, 10300, 10500, 10700, 10900};

        codingStart = 10100;
        codingEnd = 10998;
        transData = createTransExons(geneId1, transId++, POS_STRAND, exonStarts, 99, codingStart, codingEnd, true, "");
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId2, transDataList);

        // add upstream breakends
        List<BreakendGeneData> upGenes = Lists.newArrayList();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, CHR_1, 250, POS_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(1, true, CHR_1, 450, POS_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(2, true, CHR_1, 650, POS_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(3, true, CHR_1, 850, POS_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(CHR_1, 250, POS_ORIENT);
        upGenes.get(1).setPositionalData(CHR_1, 450, POS_ORIENT);
        upGenes.get(2).setPositionalData(CHR_1, 650, POS_ORIENT);
        upGenes.get(3).setPositionalData(CHR_1, 850, POS_ORIENT);

        // add downstream breakends
        List<BreakendGeneData> downGenes = Lists.newArrayList();
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, CHR_1, 10250, NEG_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(1, false, CHR_1, 10450, NEG_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(2, false, CHR_1, 10650, NEG_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(3, false, CHR_1, 10850, NEG_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(CHR_1, 10250, NEG_ORIENT);
        downGenes.get(1).setPositionalData(CHR_1, 10450, NEG_ORIENT);
        downGenes.get(2).setPositionalData(CHR_1, 10650, NEG_ORIENT);
        downGenes.get(3).setPositionalData(CHR_1, 10850, NEG_ORIENT);

        FusionParameters params = new FusionParameters();
        params.RequirePhaseMatch = true;
        params.AllowExonSkipping = true;

        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);
        fusions.forEach(x -> x.setKnownType(KNOWN_PAIR));
        tester.FusionAnalyser.getFusionFinder().findTopReportableFusion(fusions);

        assertEquals(12, fusions.size());
        final GeneFusion reportedFusion = fusions.stream().filter(x -> x.reportable()).findFirst().orElse(null);
        assertNotNull(reportedFusion);

        // the selected fusion is the longest for coding bases and without any exon skipping
        assertEquals(450, reportedFusion.upstreamTrans().gene().position());
        assertEquals(10850, reportedFusion.downstreamTrans().gene().position());
        assertEquals(0, reportedFusion.getExonsSkipped(true));
        assertEquals(0, reportedFusion.getExonsSkipped(false));
        assertTrue(reportedFusion.reportable());

        for(GeneFusion fusion: fusions)
        {
            if(fusion == reportedFusion)
                continue;

            assertTrue(!fusion.reportable());
        }
    }

    @Test
    public void testFusionPrioritisation()
    {
        // create a set of valid fusions and successively invalidate the top one to test prioritisation logic
        String geneId1 = "ENSG0001";

        BreakendGeneData upGene = createGeneAnnotation(0, true, geneId1, geneId1, POS_STRAND, CHR_1, 450, 1);

        Integer codingStart = 350;
        Integer codingEnd = 950;

        BreakendTransData upTrans1 = createTranscript(
                upGene, 1, true, 100, 1000, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING,
                2, 3, PHASE_1, 50, 300);;

        BreakendGeneData downGene = createGeneAnnotation(0, true, geneId1, geneId1, POS_STRAND, CHR_1, 450, 1);

        codingStart = 10350;
        codingEnd = 10950;

        BreakendTransData downTrans1 = createTranscript(
                downGene, 11, true, 10000, 11000, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING,
                2, 3, PHASE_1, 50, 300);

        GeneFusion fusion1 = new GeneFusion(upTrans1, downTrans1, true);
        GeneFusion fusion2 = new GeneFusion(upTrans1, downTrans1, false);

        List<GeneFusion> fusions = Lists.newArrayList(fusion1, fusion2);

        // 1. phase-matched
        GeneFusion topFusion = findTopPriorityFusion(fusions);
        assertEquals(fusion1, topFusion);

        // 2. valid chain
        fusion2 = new GeneFusion(upTrans1, downTrans1, true);

        FusionAnnotations annotations = ImmutableFusionAnnotations.builder()
                .clusterId(1).clusterCount(10).resolvedType("")
                .terminatedUp(true)
                .terminatedDown(false)
                .build();

        fusion1.setAnnotations(annotations);

        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = findTopPriorityFusion(fusions);
        assertEquals(fusion2, topFusion);

        // 3. down protein coding
        BreakendTransData downTrans2 = createTranscript(
                downGene, 12, true, 10000, 11000, codingStart, codingEnd, BIOTYPE_PROCESSED_TRANS,
                2, 3, PHASE_1, 50, 300);

        fusion1.setAnnotations(null);
        fusion2 = new GeneFusion(upTrans1, downTrans2, true);

        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = findTopPriorityFusion(fusions);
        assertEquals(fusion1, topFusion);

        // 4. exons skipped
        fusion2 = new GeneFusion(upTrans1, downTrans1, true);
        fusion1.setExonsSkipped(2, 0);
        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = findTopPriorityFusion(fusions);
        assertEquals(fusion2, topFusion);

        // 5. 3P partner canonical
        BreakendTransData downTrans3 = createTranscript(
                downGene, 13, false, 10000, 11000, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING,
                2, 3, PHASE_1, 50, 300);

        fusion1 = new GeneFusion(upTrans1, downTrans3, true);
        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = findTopPriorityFusion(fusions);
        assertEquals(fusion2, topFusion);

        // 6. 5P partner canonical
        BreakendTransData upTrans2 = createTranscript(
                upGene, 2, false, 100, 1000, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING,
                2, 3, PHASE_1, 50, 300);

        fusion1 = new GeneFusion(upTrans2, downTrans1, true);
        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = findTopPriorityFusion(fusions);
        assertEquals(fusion2, topFusion);

        // 7. 5P partner less coding bases
        BreakendTransData upTrans3 = createTranscript(
                upGene, 3, true, 100, 1000, codingStart, codingEnd, BIOTYPE_PROTEIN_CODING,
                2, 3, PHASE_1, 40, 200);

        fusion1 = new GeneFusion(upTrans3, downTrans1, true);
        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = findTopPriorityFusion(fusions);
        assertEquals(fusion2, topFusion);
    }

    @Test
    public void testChainedFusions()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        String geneName1 = "GENE1";
        String geneId1 = "ENSG0001";

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, geneName1, CHR_1, 1, 1000, 2000));

        List<TranscriptData> transDataList = Lists.newArrayList();

        String transName1 = "ENST0001";
        int transId1 = 1;

        boolean isCanonical = true;
        int transStart = 1000;
        int transEnd = 2000;
        int codingStart = 1401;
        int codingEnd = 1900;

        TranscriptData transData = new TranscriptData(transId1, transName1, geneId1, isCanonical, POS_STRAND, transStart, transEnd,
                codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        List<ExonData> exons = Lists.newArrayList();

        exons.add(new ExonData(transId1, 1000, 1100, 1, -1, -1));
        exons.add(new ExonData(transId1, 1300, 1500, 2, -1, 1));
        exons.add(new ExonData(transId1, 1600, 1700, 3, 1, 0));
        exons.add(new ExonData(transId1, 1800, 1900, 4, 0, -1));

        transData.setExons(exons);
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId1, transDataList);

        transDataList = Lists.newArrayList();

        String geneName2 = "GENE2";
        String geneId2 = "ENSG0002";

        geneList.add(createEnsemblGeneData(geneId2, geneName2, CHR_1, POS_STRAND, 10000, 12000));
        addGeneData(geneTransCache, CHR_1, geneList);

        String transName2 = "ENST0002";
        int transId2 = 2;

        transStart = 11000;
        transEnd = 12000;
        codingStart = 11049;
        codingEnd = 11980;

        transData = new TranscriptData(transId2, transName2, geneId2, isCanonical, POS_STRAND, transStart, transEnd,
                codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        exons = Lists.newArrayList();

        exons.add(new ExonData(transId1, 11000, 11100, 1, -1, 1));
        exons.add(new ExonData(transId1, 11300, 11501, 2, 1, 2));
        exons.add(new ExonData(transId1, 11600, 11699, 3, 2, 0));
        exons.add(new ExonData(transId1, 11950, 12000, 4, 0, -1));

        transData.setExons(exons);
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId2, transDataList);

        PRE_GENE_PROMOTOR_DISTANCE = 200;

        // set known fusion genes
        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache()
                .addData(new KnownFusionData(KNOWN_PAIR, geneName1, geneName2, "", ""));

        // test 1: create a chain of DELs with a single-SV fusion which link between exon 2-3 of upstream to 2-3 of downstream

        // pre-gene
        SvVarData var1 = createDel(0, CHR_1, 300,400);

        // pre-gene (just to keep the cluster big enough to not resolve into simple SVs
        SvVarData var2 = createDel(1, CHR_1, 500,600);

        // fusing del
        SvVarData var3 = createDel(2, CHR_1, 1550,11200);

        // intronic del which then runs out remainder of transcript
        SvVarData var4 = createDel(3, CHR_1, 15000,16000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        setSvGeneData(tester.AllVariants, geneTransCache, true, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        GeneFusion fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(var3.id(), fusion.upstreamTrans().gene().id());
        assertEquals(var3.id(), fusion.downstreamTrans().gene().id());
        assertTrue(validateFusionAnnotations(fusion, true, true));

        // test 2: this time a chain from the first to the last variant with the middle 2 going out to non-disruptive locations
        tester.clearClustersAndSVs();

        // upstream trans
        var1 = createDel(0, CHR_1, 1550,50000);

        var2 = createDel(1, CHR_1, 50100,51000);

        // fusing del
        var3 = createDel(2, CHR_1, 55000,56000);

        // intronic del which then runs out remainder of transcript
        var4 = createDup(3, CHR_1, 10900, 56500);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        setSvGeneData(tester.AllVariants, geneTransCache, true, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(var1.id(), fusion.upstreamTrans().gene().id());
        assertEquals(var4.id(), fusion.downstreamTrans().gene().id());
        assertTrue(validateFusionAnnotations(fusion, true, true));

        // test 4: invalid fusion, with a TI beyond the fusion ending in an exon upstream and skipping an exon downstream
        tester.clearClustersAndSVs();

        // del starting in an exon
        var1 = createDel(0, CHR_1, 1150,1250);

        // intronic dup
        var2 = createDup(1, CHR_1, 1350,1575);

        // fusion the 2 genes
        var3 = createDel(2, CHR_1, 1750,11525);

        // del skips an exon (making it disruptive to the transcript)
        var4 = createDel(3, CHR_1, 11575, 111800);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        setSvGeneData(tester.AllVariants, geneTransCache, true, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(3, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().stream().filter(x -> x.reportable()).findFirst().orElse(null);
        assertTrue(fusion != null);
        assertEquals(var3.id(), fusion.upstreamTrans().gene().id());
        assertEquals(var3.id(), fusion.downstreamTrans().gene().id());

        assertTrue(validateFusionAnnotations(fusion, false, true));

        // test 5: test a chained fusion which has invalid traversal since it goes through another gene's splice acceptor
        transDataList = Lists.newArrayList();

        String geneName3 = "GENE3";
        String geneId3 = "ENSG0003";
        String remoteChromosome = "4";

        geneList.add(createEnsemblGeneData(geneId3, geneName3, remoteChromosome, 1, 1000, 3000));
        addGeneData(geneTransCache, remoteChromosome, geneList);

        int transId3 = 3;
        int[] exonStarts = new int[]{1000, 2000, 3000};
        Integer codingStart2 = 1003;
        Integer codingEnd2 = 3098;
        TranscriptData transData3 = createTransExons(geneId3, transId3, POS_STRAND, exonStarts, 98,
                codingStart2, codingEnd2, true, "");
        transDataList.add(transData3);

        addTransExonData(geneTransCache, geneId3, transDataList);

        // test 4: invalid fusion, with a TI beyond the fusion ending in an exon upstream and skipping an exon downstream
        tester.clearClustersAndSVs();

        // inv coming in before the gene
        var1 = createBnd(1, "1", 51000, 1, "2", 1000, -1);

        // runs into start of gene
        var2 = createInv(2, CHR_1, 500, 50000, -1);

        // fuses the genes with a DEL but which goes through a few remote TIs in between including one which goes through gene 3
        var3 = createBnd(3, CHR_1, 1750, 1, "2", 10000, -1);
        var4 = createBnd(4, CHR_1, 61000, 1, "2", 11000, 1);
        SvVarData var5 = createBnd(5, CHR_1, 60000, -1, remoteChromosome, 1200, -1);
        SvVarData var6 = createBnd(6, remoteChromosome, 2500, 1, "5", 1000, -1);
        SvVarData var7 = createBnd(7, CHR_1, 11525, -1, "5", 2000, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);
        assertEquals(1, cluster.getChains().size());
        assertEquals(7, cluster.getChains().get(0).getSvCount());

        setSvGeneData(tester.AllVariants, geneTransCache, true, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(2, tester.FusionAnalyser.getFusions().size());
        assertEquals(1, tester.FusionAnalyser.getFusions().stream().filter(x -> x.validChainTraversal()).count());

        fusion = tester.FusionAnalyser.getFusions().stream().filter(x -> x.validChainTraversal()).findFirst().orElse(null);
        assertTrue(fusion.phaseMatched());
        assertEquals(var3.id(), fusion.upstreamTrans().gene().id());
        assertEquals(var5.id(), fusion.downstreamTrans().gene().id());

        assertFalse(validateFusionAnnotations(fusion, false, false));

        // repeat again but this time allowing for invalid traversal in a known pair
        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(2, tester.FusionAnalyser.getFusions().size());

        final int varId3 = var3.id();
        final int varId7 = var7.id();

        fusion = tester.FusionAnalyser.getFusions().stream()
                .filter(x -> x.upstreamTrans().gene().id() == varId3 && x.downstreamTrans().gene().id() == varId7).findFirst().orElse(null);

        assertNotNull(fusion);
        assertEquals(KNOWN_PAIR, fusion.knownType());
        assertFalse(fusion.validChainTraversal());


        // test 5: non-disruptive chain, doesn't register a fusion
        PRE_GENE_PROMOTOR_DISTANCE = 30000;

        tester.clearClustersAndSVs();

        // inv coming in before the gene
        var1 = createBnd(1, "1", 1510, 1, "2", 1000, -1);
        var2 = createBnd(2, "1", 1580, -1, "2", 1100, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);
        assertEquals(1, cluster.getChains().size());

        setSvGeneData(tester.AllVariants, geneTransCache, true, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertTrue(tester.FusionAnalyser.getFusions().isEmpty());
        assertEquals(1, tester.FusionAnalyser.getInvalidFusions().size());
        assertTrue(tester.FusionAnalyser.getInvalidFusions().keySet().iterator().next().nonDisruptiveChain());


        // test 6: a chained fusion with non-disrputive SVs at both ends, means the fusion is not chain terminated
        tester.clearClustersAndSVs();

        // fully intronic, non-disruptive del
        var1 = createDel(0, CHR_1, 1150,1250);

        // intronic dup
        var2 = createDup(1, CHR_1, 1525,1575);

        // fuses the 2 genes
        var3 = createDel(2, CHR_1, 1750,11525);

        // fully intronic, non-disruptive TI
        var4 = createBnd(3, CHR_1, 11750, 1, CHR_2, 50000, -1);
        var5 = createBnd(4, CHR_1, 11850, -1, CHR_2, 50100, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());
        SvChain chain = cluster.getChains().get(0);
        assertEquals(4, chain.getLinkCount());

        setSvGeneData(tester.AllVariants, geneTransCache, true, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getUniqueFusions().size());

        fusion = tester.FusionAnalyser.getUniqueFusions().stream().filter(x -> x.reportable()).findFirst().orElse(null);
        assertTrue(fusion != null);
        assertEquals(var3.id(), fusion.upstreamTrans().gene().id());
        assertEquals(var3.id(), fusion.downstreamTrans().gene().id());

        assertTrue(validateFusionAnnotations(fusion, true, true));
    }

    private static boolean validateFusionAnnotations(final GeneFusion fusion, boolean validEnds, boolean validTraversal)
    {
        final FusionAnnotations annotations = fusion.getAnnotations();
        if(annotations == null)
            return false;

        final FusionChainInfo chainInfo = annotations.chainInfo();

        if(chainInfo == null)
            return false;

        if(validTraversal != chainInfo.validTraversal())
            return false;

        // test the exons disrupted and terminated fields
        boolean validUp = !annotations.terminatedUp();
        boolean validDown = !annotations.terminatedDown();

        if(validEnds)
            return !fusion.isTerminated();
        else
            return !validUp || !validDown;
    }

}
