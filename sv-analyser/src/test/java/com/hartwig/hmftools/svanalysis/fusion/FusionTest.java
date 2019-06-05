package com.hartwig.hmftools.svanalysis.fusion;

import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_KNOWN;
import static com.hartwig.hmftools.svanalysis.analyser.com.hartwig.hmftools.svanalysis.gene.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.svanalysis.analyser.com.hartwig.hmftools.svanalysis.gene.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.svanalysis.analyser.com.hartwig.hmftools.svanalysis.gene.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.svanalysis.fusion.FusionDisruptionAnalyser.FCI_TRAV_ASSEMBLY;
import static com.hartwig.hmftools.svanalysis.fusion.FusionDisruptionAnalyser.FCI_VALID_TRAVERSAL;
import static com.hartwig.hmftools.svanalysis.fusion.FusionDisruptionAnalyser.isDisrupted;
import static com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection.PRE_GENE_PROMOTOR_DISTANCE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.svanalysis.analyser.SvTestHelper;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection;

import org.junit.Test;

public class FusionTest
{
    @Test
    public void testReportableFusionComparison()
    {
        // test the selection and prioritisation logic from a collection of valid fusions
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        SvGeneTranscriptCollection geneTransCache = new SvGeneTranscriptCollection();

        tester.initialiseFusions(geneTransCache);

        PRE_GENE_PROMOTOR_DISTANCE = 100;

        // first a gene on the forward strand
        String geneName = "GENE1";
        String geneId1 = "ENSG0001";
        String chromosome = "1";
        byte strand = 1;

        // GeneAnnotation geneUp = createGeneAnnotation(0, true, geneName, geneId1, strand, chromosome, 0, 1);

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, geneName, chromosome, 1, 100, 1000));

        List<TranscriptExonData> transExonList = Lists.newArrayList();

        String transName = "ENST0001";
        int transId = 1;

        long transStart = 100;
        long transEnd = 1000;
        long codingStart = 350;
        long codingEnd = 950;

        transExonList.add(new TranscriptExonData(geneId1, transName, transId, true, strand, transStart, transEnd,
                100, 200, 1, -1, -1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId1, transName, transId, true, strand, transStart, transEnd,
                300, 400, 2, -1, 1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId1, transName, transId, true, strand, transStart, transEnd,
                500, 600, 3, 1, 2, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId1, transName, transId, true, strand, transStart, transEnd,
                700, 800, 4, 2, 0, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId1, transName, transId, true, strand, transStart, transEnd,
                900, 1000, 5, 0, -1, codingStart, codingEnd, ""));

        addTransExonData(geneTransCache, geneId1, transExonList);

        geneName = "GENE2";
        String geneId2 = "ENSG0002";

        // GeneAnnotation geneDown = createGeneAnnotation(0, true, geneName, geneId1, strand, chromosome, 0, -1);

        geneList.add(createEnsemblGeneData(geneId2, geneName, chromosome, 1, 10000, 12000));

        addGeneData(geneTransCache, chromosome, geneList);

        String transName2 = "ENST0002";
        int transId2 = 2;

        transStart = 10100;
        transEnd = 11000;
        codingStart = 10150;
        codingEnd = 10750;

        transExonList = Lists.newArrayList();

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, true, strand, transStart, transEnd,
                10100, 10200, 1, -1, 1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, true, strand, transStart, transEnd,
                10300, 10400, 2, 1, 2, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, true, strand, transStart, transEnd,
                10500, 10600, 3, 2, 0, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, true, strand, transStart, transEnd,
                10700, 10800, 4, 0, -1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, true, strand, transStart, transEnd,
                10900, 11000, 5, -1, -1, codingStart, codingEnd, ""));

        addTransExonData(geneTransCache, geneId2, transExonList);

        byte posOrient = 1;
        byte negOrient = -1;

        // add upstream breakends
        List<GeneAnnotation> upGenes = Lists.newArrayList();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome, 250, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(1, true, chromosome, 450, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(2, true, chromosome, 650, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(3, true, chromosome, 850, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome, 250, posOrient);
        upGenes.get(1).setPositionalData(chromosome, 450, posOrient);
        upGenes.get(2).setPositionalData(chromosome, 650, posOrient);
        upGenes.get(3).setPositionalData(chromosome, 850, posOrient);

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

        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, true, true, null, false);
        fusions.forEach(x -> x.setKnownFusionType(REPORTABLE_TYPE_KNOWN));
        tester.FusionAnalyser.getFusionFinder().setReportableGeneFusions(fusions);

        assertEquals(6, fusions.size());
        final GeneFusion fusion = fusions.get(0);

        // the selected fusion is the longest for coding bases and without any exon skipping
        assertEquals(450, fusion.upstreamTrans().parent().position());
        assertEquals(10250, fusion.downstreamTrans().parent().position());
        assertEquals(0, fusion.getExonsSkipped(true));
        assertEquals(0, fusion.getExonsSkipped(false));
        assertTrue(fusion.reportable());

        for(int i = 1; i < fusions.size(); ++i)
        {
            assertTrue(!fusions.get(i).reportable());
        }
    }

    @Test
    public void testChainedFusions()
    {
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        SvGeneTranscriptCollection geneTransCache = new SvGeneTranscriptCollection();

        tester.initialiseFusions(geneTransCache);

        String geneName1 = "GENE1";
        String geneId1 = "ENSG0001";
        String chromosome = "1";

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, geneName1, chromosome, 1, 1000, 2000));

        List<TranscriptExonData> transExonList = Lists.newArrayList();

        String transName1 = "ENST0001";
        int transId1 = 1;

        byte strand = (byte)1;
        boolean isCanonical = true;
        long transStart = 1000;
        long transEnd = 2000;
        long codingStart = 1400;
        long codingEnd = 1900;
        transExonList.add(new TranscriptExonData(geneId1, transName1, transId1, isCanonical, strand, transStart, transEnd,
                1000, 1100, 1, -1, -1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId1, transName1, transId1, isCanonical, strand, transStart, transEnd,
                1300, 1500, 2, -1, 1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId1, transName1, transId1, isCanonical, strand, transStart, transEnd,
                1600, 1700, 3, 1, 2, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId1, transName1, transId1, isCanonical, strand, transStart, transEnd,
                1800, 1900, 4, 2, -1, codingStart, codingEnd, ""));

        addTransExonData(geneTransCache, geneId1, transExonList);

        String geneName2 = "GENE2";
        String geneId2 = "ENSG0002";

        geneList.add(createEnsemblGeneData(geneId2, geneName2, chromosome, 1, 10000, 12000));
        addGeneData(geneTransCache, chromosome, geneList);

        String transName2 = "ENST0002";
        int transId2 = 2;

        transStart = 11000;
        transEnd = 12000;
        codingStart = 11050;
        codingEnd = 11980;

        transExonList = Lists.newArrayList();

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, isCanonical, strand, transStart, transEnd,
                11000, 11100, 1, -1, 1, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, isCanonical, strand, transStart, transEnd,
                11300, 11500, 2, 1, 2, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, isCanonical, strand, transStart, transEnd,
                11600, 11700, 3, 2, 0, codingStart, codingEnd, ""));

        transExonList.add(new TranscriptExonData(geneId2, transName2, transId2, isCanonical, strand, transStart, transEnd,
                11950, 12000, 4, 2, -1, codingStart, codingEnd, ""));

        addTransExonData(geneTransCache, geneId2, transExonList);

        PRE_GENE_PROMOTOR_DISTANCE = 200;

        // test 1: create a chain of DELs with a single-SV fusion which link between exon 2-3 of upstream to 2-3 of downstream

        // pre-gene
        SvVarData var1 = createDel("0", chromosome, 300,400);

        // pre-gene (just to keep the cluster big enough to not resolve into simple SVs
        SvVarData var2 = createDel("1", chromosome, 500,600);

        // fusing del
        SvVarData var3 = createDel("2", chromosome, 1550,11200);

        // intronic del which then runs out remainder of transcript
        SvVarData var4 = createDel("3", chromosome, 15000,16000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        tester.FusionAnalyser.setSvGeneData(tester.AllVariants, geneTransCache, true, false, true);
        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.ClusteringMethods.getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        GeneFusion fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(var3.dbId(), fusion.upstreamTrans().parent().id());
        assertEquals(var3.dbId(), fusion.downstreamTrans().parent().id());
        assertTrue(validateFusionAnnotations(fusion, true, true));

        // test 2: this time a chain from the first to the last variant with the middle 2 going out to non-disruptive locations
        tester.clearClustersAndSVs();

        // upstream trans
        var1 = createDel("0", chromosome, 1550,50000);

        var2 = createDel("1", chromosome, 50100,51000);

        // fusing del
        var3 = createDel("2", chromosome, 55000,56000);

        // intronic del which then runs out remainder of transcript
        var4 = createDup("3", chromosome, 10900, 56500);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        tester.FusionAnalyser.setSvGeneData(tester.AllVariants, geneTransCache, true, false, true);
        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.ClusteringMethods.getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(var1.dbId(), fusion.upstreamTrans().parent().id());
        assertEquals(var4.dbId(), fusion.downstreamTrans().parent().id());
        assertTrue(validateFusionAnnotations(fusion, true, true));

        // test 4: invalid fusion, with a TI beyond the fusion ending in an exon upstream and skipping an exon downstream
        tester.clearClustersAndSVs();

        // del starting in an exon
        var1 = createDel("0", chromosome, 1150,1250);

        // intronic dup
        var2 = createDup("1", chromosome, 1350,1575);

        // fusion the 2 genes
        var3 = createDel("2", chromosome, 1750,11525);

        // del skips an exon
        var4 = createDel("3", chromosome, 11575, 111800);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        tester.FusionAnalyser.setSvGeneData(tester.AllVariants, geneTransCache, true, false, true);
        tester.FusionAnalyser.run(tester.SampleId, tester.AllVariants, null,
                tester.getClusters(), tester.ClusteringMethods.getChrBreakendMap());

        assertEquals(1, tester.FusionAnalyser.getFusions().size());

        fusion = tester.FusionAnalyser.getFusions().get(0);
        assertEquals(var3.dbId(), fusion.upstreamTrans().parent().id());
        assertEquals(var3.dbId(), fusion.downstreamTrans().parent().id());

        assertTrue(validateFusionAnnotations(fusion, false, true));
    }

    private static boolean validateFusionAnnotations(final GeneFusion fusion, boolean validEnds, boolean validTraversal)
    {
        String[] fields = fusion.getAnnotations().split(",");

        // PhaseMatched,ClusterId,ClusterCount,ResolvedType,OverlapUp,OverlapDown,ChainInfo
        if(fields.length != 7)
            return false;

        String[] chainInfo = fields[6].split(";");

        int fieldLength = FCI_TRAV_ASSEMBLY+1;
        if(chainInfo.length != fieldLength)
            return false;

        if(validTraversal != (chainInfo[FCI_VALID_TRAVERSAL].equals("true")))
            return false;

        String disruptionsUp = fields[4];
        String disruptionsDown = fields[5];

        // test the exons disrupted and terminated fields
        boolean validUp = !isDisrupted(disruptionsUp);
        boolean validDown = !isDisrupted(disruptionsDown);

        if(validEnds)
            return validUp && validDown;
        else
            return !validUp || !validDown;
    }

}
