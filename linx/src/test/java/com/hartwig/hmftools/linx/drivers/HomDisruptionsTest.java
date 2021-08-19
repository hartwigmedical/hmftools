package com.hartwig.hmftools.linx.drivers;

import static com.hartwig.hmftools.linx.analysis.VariantPrep.setSvGeneData;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.HOM_DEL_DISRUPTION;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.HOM_DUP_DISRUPTION;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.createDriverGene;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;

import static org.junit.Assert.assertEquals;

import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.SampleAnalyser;
import com.hartwig.hmftools.linx.analysis.VariantPrep;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class HomDisruptionsTest
{
    @Test
    public void testSimpleDelDisruption()
    {
        LinxTester tester = new LinxTester();

        tester.setNonClusterAllelePloidies(1, 0);

        EnsemblDataCache geneTransCache = createGeneDataCache();
        tester.initialiseFusions(geneTransCache);

        // must be a known TSG driver
        String geneName = "TP53";
        String geneId = "ENSG00000141510";
        String chromosome = "17";
        byte strand = 1;

        int transStart = 10000;
        int transEnd = 130000;

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, transStart, transEnd));

        tester.Config.DriverGenes.add(createDriverGene(geneName));

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;
        int[] exonStarts = new int[] { 10000, 40000, 70000, 100000, 130000};

        TranscriptData transData = createTransExons(geneId, transId++, strand, exonStarts, 100,
                10050, 130050, true, "");
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId, transDataList);
        addGeneData(geneTransCache, chromosome, geneList);

        tester.initialiseDriverGeneAnnotator(geneTransCache);
        DriverGeneAnnotator driverAnnotator = tester.DriverAnnotator;
        driverAnnotator.setSamplePurityData(2, false);

        // one pair of DELs form a DB but the DELs themselves are not disruptive and since there's no way of knowing if they're phased,
        // they aren't considered to make a hom-disruption
        SvVarData del1 = createDel(tester.nextVarId(), chromosome, 12000, 24000);
        SvVarData del2 = createDel(tester.nextVarId(), chromosome, 18000, 30000);

        SvVarData del3 = createDel(tester.nextVarId(), chromosome, 42000, 85000);
        SvVarData del4 = createDel(tester.nextVarId(), chromosome, 75000, 120000);

        tester.AllVariants.add(del1);
        tester.AllVariants.add(del2);
        tester.AllVariants.add(del3);
        tester.AllVariants.add(del4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(4, tester.Analyser.getClusters().size());

        setSvGeneData(tester.AllVariants, geneTransCache, false, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, false);

        driverAnnotator.annotateSVs(tester.SampleId, tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, driverAnnotator.getDriverGeneDataList().size());
        DriverGeneData dgData = driverAnnotator.getDriverGeneDataList().get(0);

        assertEquals(1, dgData.getEvents().size());
        assertEquals(transData, dgData.TransData);
        assertEquals(ResolvedType.DEL, dgData.getEvents().get(0).getCluster().getResolvedType());
        assertTrue(dgData.getEvents().get(0).getCluster().getSVs().contains(del4));
        assertEquals(HOM_DEL_DISRUPTION, dgData.getEvents().get(0).Type);
    }

    @Test
    public void testChainedDelDisruption()
    {
        LinxTester tester = new LinxTester();

        tester.setNonClusterAllelePloidies(0.2, 0);

        EnsemblDataCache geneTransCache = createGeneDataCache();
        tester.initialiseFusions(geneTransCache);

        // must be a known TSG driver
        String geneName = "TP53";
        String geneId = "ENSG00000141510";
        String chromosome = "17";
        byte strand = 1;

        int transStart = 10000;
        int transEnd = 55000;

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, transStart, transEnd));

        tester.Config.DriverGenes.add(createDriverGene(geneName));

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;
        int[] exonStarts = new int[] { 10000, 20000, 30000, 40000, 50000};
        // int[] exonPhases = new int[] {-1, 0, 0, 0, -1};

        int codingStart = 10002;
        int codingEnd = 50098;
        TranscriptData transData = createTransExons(geneId, transId++, strand, exonStarts, 100, codingStart, codingEnd, true, "");
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId, transDataList);
        addGeneData(geneTransCache, chromosome, geneList);

        tester.initialiseDriverGeneAnnotator(geneTransCache);
        DriverGeneAnnotator driverAnnotator = tester.DriverAnnotator;
        driverAnnotator.setSamplePurityData(2, false);

        // make a chain of SVs which contain one or more deletion bridges affecting the gene
        SvVarData var1 = createInv(tester.nextVarId(), chromosome, 15000, 25000, 1);
        SvVarData var2 = createDup(tester.nextVarId(), chromosome, 15100, 35000);
        SvVarData var3 = createInv(tester.nextVarId(), chromosome, 25100, 35100, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        setSvGeneData(tester.AllVariants, geneTransCache, false, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, false);

        driverAnnotator.annotateSVs(tester.SampleId, tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, driverAnnotator.getDriverGeneDataList().size());
        DriverGeneData dgData = driverAnnotator.getDriverGeneDataList().get(0);

        assertEquals(1, dgData.getEvents().size());
        assertEquals(cluster, dgData.getEvents().get(0).getCluster());
        assertEquals(HOM_DEL_DISRUPTION, dgData.getEvents().get(0).Type);

        tester.clearClustersAndSVs();

        // add an unrelated DEL to set the correct copy number running out to telomere
        SvVarData del = createDel(tester.nextVarId(), chromosome, 400, 500);

        var1 = createInv(tester.nextVarId(), chromosome, 15000, 25000, 1);
        var2 = createInv(tester.nextVarId(), chromosome, 14990, 24990, -1);

        tester.AllVariants.add(del);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        setSvGeneData(tester.AllVariants, geneTransCache, false, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, false);

        driverAnnotator.clearResults(true);
        driverAnnotator.annotateSVs(tester.SampleId, tester.Analyser.getState().getChrBreakendMap());

        assertEquals(2, tester.Analyser.getClusters().size());
        cluster = tester.findClusterWithSVs(Lists.newArrayList(var1, var2));

        assertEquals(1, driverAnnotator.getDriverGeneDataList().size());
        dgData = driverAnnotator.getDriverGeneDataList().get(0);

        assertEquals(1, dgData.getEvents().size());
        assertEquals(cluster, dgData.getEvents().get(0).getCluster());
        assertEquals(HOM_DEL_DISRUPTION, dgData.getEvents().get(0).Type);
    }

    @Test
    public void testDupDisruption()
    {
        LinxTester tester = new LinxTester();

        tester.setNonClusterAllelePloidies(0.2, 0);

        EnsemblDataCache geneTransCache = createGeneDataCache();
        tester.initialiseFusions(geneTransCache);

        // must be a known TSG driver
        String geneName = "TP53";
        String geneId = "ENSG00000141510";
        String chromosome = "17";
        byte strand = 1;

        int transStart = 10000;
        int transEnd = 55000;

        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, transStart, transEnd));

        tester.Config.DriverGenes.add(createDriverGene(geneName));

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;
        int[] exonStarts = new int[] { 10000, 20000, 30000, 40000, 50000};
        // int[] exonPhases = new int[] {-1, 0, 0, 0, -1};

        int codingStart = 10002;
        int codingEnd = 50098;
        TranscriptData transData = createTransExons(geneId, transId++, strand, exonStarts, 100, codingStart, codingEnd, true, "");
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId, transDataList);
        addGeneData(geneTransCache, chromosome, geneList);

        tester.initialiseDriverGeneAnnotator(geneTransCache);
        DriverGeneAnnotator driverAnnotator = tester.DriverAnnotator;
        driverAnnotator.setSamplePurityData(2, false);

        // the DUP must be disruptive, ie duplicating at least 1 exon
        SvVarData var1 = createDup(tester.nextVarId(), chromosome, 15000, 35000);

        tester.AllVariants.add(var1);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        setSvGeneData(tester.AllVariants, geneTransCache, false, false);
        tester.FusionAnalyser.annotateTranscripts(tester.AllVariants, false);

        driverAnnotator.annotateSVs(tester.SampleId, tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, driverAnnotator.getDriverGeneDataList().size());
        DriverGeneData dgData = driverAnnotator.getDriverGeneDataList().get(0);

        assertEquals(1, dgData.getEvents().size());
        assertEquals(cluster, dgData.getEvents().get(0).getCluster());
        assertEquals(HOM_DUP_DISRUPTION, dgData.getEvents().get(0).Type);
    }

}
