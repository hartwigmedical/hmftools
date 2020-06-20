package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.generateExonStarts;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROCESSED_TRANS;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneAnnotation;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.determineReportableFusion;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createTransExons;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.linx.analysis.SvSampleAnalyser;
import com.hartwig.hmftools.linx.utils.LinxTester;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

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
        String chromosome = "1";
        byte strand = 1;

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, geneName, chromosome, strand, 100, 1000));

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;

        int[] exonStarts = new int[]{100, 300, 500, 700, 900};
        int[] exonPhases = new int[]{-1, 1, 2, 0, -1};

        TranscriptData transData = createTransExons(geneId1, transId++, strand, exonStarts, exonPhases, 100, true);
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId1, transDataList);

        geneName = "GENE2";
        String geneId2 = "ENSG0002";

        // GeneAnnotation geneDown = createGeneAnnotation(0, true, geneName, geneId1, strand, chromosome, 0, -1);

        geneList.add(createEnsemblGeneData(geneId2, geneName, chromosome, 1, 10000, 12000));

        addGeneData(geneTransCache, chromosome, geneList);

        transDataList = Lists.newArrayList();

        exonStarts = new int[]{10100, 10300, 10500, 10700, 10900};
        exonPhases = new int[]{1, 2, 0, -1, -1};

        transData = createTransExons(geneId1, transId++, strand, exonStarts, exonPhases, 100, true);
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId2, transDataList);

        byte POS_ORIENT = 1;
        byte negOrient = -1;

        // add upstream breakends
        List<GeneAnnotation> upGenes = Lists.newArrayList();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome, 250, POS_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(1, true, chromosome, 450, POS_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(2, true, chromosome, 650, POS_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(3, true, chromosome, 850, POS_ORIENT, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome, 250, POS_ORIENT);
        upGenes.get(1).setPositionalData(chromosome, 450, POS_ORIENT);
        upGenes.get(2).setPositionalData(chromosome, 650, POS_ORIENT);
        upGenes.get(3).setPositionalData(chromosome, 850, POS_ORIENT);

        // add downstream breakends
        List<GeneAnnotation> downGenes = Lists.newArrayList();
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome, 10250, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(1, false, chromosome, 10450, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(2, false, chromosome, 10650, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(3, false, chromosome, 10850, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome, 10250, negOrient);
        downGenes.get(1).setPositionalData(chromosome, 10450, negOrient);
        downGenes.get(2).setPositionalData(chromosome, 10650, negOrient);
        downGenes.get(3).setPositionalData(chromosome, 10850, negOrient);

        FusionParameters params = new FusionParameters();
        params.RequirePhaseMatch = true;
        params.AllowExonSkipping = true;

        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params, false);
        fusions.forEach(x -> x.setKnownType(KNOWN_PAIR));
        tester.FusionAnalyser.getFusionFinder().setReportableGeneFusions(fusions);

        assertEquals(6, fusions.size());
        final GeneFusion fusion = fusions.get(0);

        // the selected fusion is the longest for coding bases and without any exon skipping
        assertEquals(450, fusion.upstreamTrans().gene().position());
        assertEquals(10250, fusion.downstreamTrans().gene().position());
        assertEquals(0, fusion.getExonsSkipped(true));
        assertEquals(0, fusion.getExonsSkipped(false));
        assertTrue(fusion.reportable());

        for(int i = 1; i < fusions.size(); ++i)
        {
            assertTrue(!fusions.get(i).reportable());
        }
    }

    @Test
    public void testFusionPrioritisation()
    {
        // create a set of valid fusions and successively invalidate the top one to test prioritisation logic
        String chromosome = "1";

        String geneId1 = "ENSG0001";

        GeneAnnotation upGene = createGeneAnnotation(0, true, geneId1, geneId1, 1, chromosome, 450, 1);

        Integer codingStart = new Integer(350);
        Integer codingEnd = new Integer(950);

        Transcript upTrans1 = new Transcript(upGene, 1, "TRANS01", 2, 1, 3, 1,
                50, 300,5, true, 100, 1000, codingStart, codingEnd);
        upTrans1.setBioType(BIOTYPE_PROTEIN_CODING);

        GeneAnnotation downGene = createGeneAnnotation(0, true, geneId1, geneId1, 1, chromosome, 450, 1);

        codingStart = new Integer(10350);
        codingEnd = new Integer(10950);

        Transcript downTrans1 = new Transcript(downGene, 11, "TRANS11", 2, 1, 3, 1,
                50, 300,5, true, 10000, 11000, codingStart, codingEnd);
        downTrans1.setBioType(BIOTYPE_PROTEIN_CODING);

        GeneFusion fusion1 = new GeneFusion(upTrans1, downTrans1, true);
        GeneFusion fusion2 = new GeneFusion(upTrans1, downTrans1, false);

        List<GeneFusion> fusions = Lists.newArrayList(fusion1, fusion2);

        // 1. phase-matched
        GeneFusion topFusion = determineReportableFusion(fusions, false);
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
        topFusion = determineReportableFusion(fusions, false);
        assertEquals(fusion2, topFusion);

        // 3. down protein coding
        Transcript downTrans2 = new Transcript(downGene, 12, "TRANS12", 2, 1, 3, 1,
                50, 300,5, true, 10000, 11000, codingStart, codingEnd);
        downTrans2.setBioType(BIOTYPE_PROCESSED_TRANS);

        fusion1.setAnnotations(null);
        fusion2 = new GeneFusion(upTrans1, downTrans2, true);

        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = determineReportableFusion(fusions, false);
        assertEquals(fusion1, topFusion);

        // 4. exons skipped
        fusion2 = new GeneFusion(upTrans1, downTrans1, true);
        fusion1.setExonsSkipped(2, 0);
        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = determineReportableFusion(fusions, false);
        assertEquals(fusion2, topFusion);

        // 5. 3P partner canonical
        Transcript downTrans3 = new Transcript(downGene, 13, "TRANS13", 2, 1, 3, 1,
                50, 300,5, false, 10000, 11000, codingStart, codingEnd);
        downTrans3.setBioType(BIOTYPE_PROTEIN_CODING);

        fusion1 = new GeneFusion(upTrans1, downTrans3, true);
        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = determineReportableFusion(fusions, false);
        assertEquals(fusion2, topFusion);

        // 6. 5P partner canonical
        Transcript upTrans2 = new Transcript(upGene, 2, "TRANS02", 2, 1, 3, 1,
                50, 300,5, false, 100, 1000, codingStart, codingEnd);
        upTrans2.setBioType(BIOTYPE_PROTEIN_CODING);

        fusion1 = new GeneFusion(upTrans2, downTrans1, true);
        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = determineReportableFusion(fusions, false);
        assertEquals(fusion2, topFusion);

        // 7. 5P partner less coding bases
        Transcript upTrans3 = new Transcript(upGene, 3, "TRANS03", 2, 1, 3, 1,
                40, 200,5, true, 100, 1000, codingStart, codingEnd);
        upTrans3.setBioType(BIOTYPE_PROTEIN_CODING);

        fusion1 = new GeneFusion(upTrans3, downTrans1, true);
        fusions = Lists.newArrayList(fusion1, fusion2);
        topFusion = determineReportableFusion(fusions, false);
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
        String chromosome = "1";

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, geneName1, chromosome, 1, 1000, 2000));

        List<TranscriptData> transDataList = Lists.newArrayList();

        String transName1 = "ENST0001";
        int transId1 = 1;

        byte strand = POS_STRAND;
        boolean isCanonical = true;
        int transStart = 1000;
        int transEnd = 2000;
        int codingStart = 1400;
        int codingEnd = 1900;

        TranscriptData transData = new TranscriptData(transId1, transName1, geneId1, isCanonical, strand, transStart, transEnd,
                codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        List<ExonData> exons = Lists.newArrayList();

        exons.add(new ExonData(transId1, 1000, 1100, 1, -1, -1));
        exons.add(new ExonData(transId1, 1300, 1500, 2, -1, 1));
        exons.add(new ExonData(transId1, 1600, 1700, 3, 1, 2));
        exons.add(new ExonData(transId1, 1800, 1900, 4, 2, -1));

        transData.setExons(exons);
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId1, transDataList);

        transDataList = Lists.newArrayList();

        String geneName2 = "GENE2";
        String geneId2 = "ENSG0002";

        geneList.add(createEnsemblGeneData(geneId2, geneName2, chromosome, strand, 10000, 12000));
        addGeneData(geneTransCache, chromosome, geneList);

        String transName2 = "ENST0002";
        int transId2 = 2;

        transStart = 11000;
        transEnd = 12000;
        codingStart = 11050;
        codingEnd = 11980;

        transData = new TranscriptData(transId2, transName2, geneId2, isCanonical, strand, transStart, transEnd,
                codingStart, codingEnd, BIOTYPE_PROTEIN_CODING);

        exons = Lists.newArrayList();

        exons.add(new ExonData(transId1, 11000, 11100, 1, -1, 1));
        exons.add(new ExonData(transId1, 11300, 11500, 2, 1, 2));
        exons.add(new ExonData(transId1, 11600, 11700, 3, 2, 0));
        exons.add(new ExonData(transId1, 11950, 12000, 4, 2, -1));

        transData.setExons(exons);
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId2, transDataList);

        PRE_GENE_PROMOTOR_DISTANCE = 200;

        // set known fusion genes
        tester.FusionAnalyser.getFusionFinder().getKnownFusionCache()
                .addData(new KnownFusionData(KNOWN_PAIR, geneName1, geneName2, "", "", ""));

        // test 1: create a chain of DELs with a single-SV fusion which link between exon 2-3 of upstream to 2-3 of downstream

        // pre-gene
        SvVarData var1 = createDel(0, chromosome, 300,400);

        // pre-gene (just to keep the cluster big enough to not resolve into simple SVs
        SvVarData var2 = createDel(1, chromosome, 500,600);

        // fusing del
        SvVarData var3 = createDel(2, chromosome, 1550,11200);

        // intronic del which then runs out remainder of transcript
        SvVarData var4 = createDel(3, chromosome, 15000,16000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        SvSampleAnalyser.setSvGeneData(tester.AllVariants, geneTransCache, true, false);
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
        var1 = createDel(0, chromosome, 1550,50000);

        var2 = createDel(1, chromosome, 50100,51000);

        // fusing del
        var3 = createDel(2, chromosome, 55000,56000);

        // intronic del which then runs out remainder of transcript
        var4 = createDup(3, chromosome, 10900, 56500);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        SvSampleAnalyser.setSvGeneData(tester.AllVariants, geneTransCache, true, false);
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
        var1 = createDel(0, chromosome, 1150,1250);

        // intronic dup
        var2 = createDup(1, chromosome, 1350,1575);

        // fusion the 2 genes
        var3 = createDel(2, chromosome, 1750,11525);

        // del skips an exon
        var4 = createDel(3, chromosome, 11575, 111800);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        SvSampleAnalyser.setSvGeneData(tester.AllVariants, geneTransCache, true, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().get(0);
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
        int[] exonPhases = new int[]{0, 0, 0,};
        TranscriptData transData3 = createTransExons(geneId3, transId3, strand, exonStarts, exonPhases, 100);
        transDataList.add(transData3);

        addTransExonData(geneTransCache, geneId3, transDataList);

        // test 4: invalid fusion, with a TI beyond the fusion ending in an exon upstream and skipping an exon downstream
        tester.clearClustersAndSVs();

        // inv coming in before the gene
        var1 = createBnd(1, "1", 51000, 1, "2", 1000, -1);

        // runs into start of gene
        var2 = createInv(2, chromosome, 500, 50000, -1);

        // fuses the genes with a DEL but which goes through a few remote TIs in between including one which goes through gene 3
        var3 = createBnd(3, chromosome, 1750, 1, "2", 10000, -1);
        var4 = createBnd(4, chromosome, 61000, 1, "2", 11000, 1);
        SvVarData var5 = createBnd(5, chromosome, 60000, -1, remoteChromosome, 1200, -1);
        SvVarData var6 = createBnd(6, remoteChromosome, 2500, 1, "5", 1000, -1);
        SvVarData var7 = createBnd(7, chromosome, 11525, -1, "5", 2000, 1);
        // var3 = createDel(2, chromosome, 1750,11525);

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

        SvSampleAnalyser.setSvGeneData(tester.AllVariants, geneTransCache, true, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, true);

        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.Analyser.getState().getChrBreakendMap());

        // invalid fusions are no longer cached
        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().get(0);
        assertTrue(!fusion.phaseMatched());
        assertEquals(var3.id(), fusion.upstreamTrans().gene().id());
        assertEquals(var5.id(), fusion.downstreamTrans().gene().id());

        assertFalse(validateFusionAnnotations(fusion, false, false));
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
