package com.hartwig.hmftools.linx.drivers;

import static com.hartwig.hmftools.common.drivercatalog.DriverType.AMP;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.linx.utils.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDriver;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createGeneCopyNumber;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createTestSv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class AmpDriversTest
{
    @Test
    public void testAmpDrivers()
    {
        LinxTester tester = new LinxTester();

        tester.setNonClusterAllelePloidies(1, 0);

        EnsemblDataCache geneTransCache = createGeneDataCache();

        DriverGeneAnnotator driverAnnotator = new DriverGeneAnnotator(null, geneTransCache, tester.Config, tester.CnDataLoader);

        driverAnnotator.setSamplePurityData(2, false);

        String geneName = GENE_NAME_1;
        String geneId = GENE_ID_1;
        String chromosome = CHR_1;
        byte strand = 1;

        int transStart = 100000;
        int transEnd = 120000;

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId, geneName, chromosome, strand, transStart, transEnd));

        List<TranscriptData> transDataList = Lists.newArrayList();

        int transId = 1;

        TranscriptData transData = new TranscriptData(transId, String.format("TRAN%04d", transId), geneId, true, strand, transStart, transEnd,
            null, null, "");
        transDataList.add(transData);

        addTransExonData(geneTransCache, geneId, transDataList);
        addGeneData(geneTransCache, chromosome, geneList);

        double geneMinCopyNumber = 12;

        final DriverCatalog driver = createDriver(geneName, chromosome, AMP, DriverCategory.ONCO, false, geneMinCopyNumber);
        final GeneCopyNumber geneCopyNumber = createGeneCopyNumber(geneName, chromosome, geneMinCopyNumber, transStart, transEnd);

        driverAnnotator.addDriverGene(driver, geneCopyNumber);

        // preceding DEL to set arm copy number
        SvVarData var1 = createDel(tester.nextVarId(), chromosome, 500, 600);

        // single DUP around the gene
        SvVarData var2 = createTestSv(tester.nextVarId(), chromosome, chromosome, transStart - 1000,transEnd + 1000,
                -1, 1, DUP,  10);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        driverAnnotator.clearResults(false);
        driverAnnotator.annotateSVs(tester.SampleId, tester.Analyser.getState().getChrBreakendMap());

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(var2));
        assertTrue(cluster != null);

        assertEquals(1, driverAnnotator.getDriverGeneDataList().size());
        DriverGeneData dgData = driverAnnotator.getDriverGeneDataList().get(0);

        assertEquals(1, dgData.getEvents().size());
        assertEquals(cluster, dgData.getEvents().get(0).getCluster());

        tester.clearClustersAndSVs();

        // on lower side: BND identified, DEL ignored, DUP ignored, INV identified
        var1 = createDel(tester.nextVarId(), chromosome, 500, 600);

        var2 = createTestSv(tester.nextVarId(), chromosome, "2", 10000,1000, -1, 1, BND,  5);

        // DEL and DUP will be ignored since are net neutral
        SvVarData var3 = createDel(tester.nextVarId(), chromosome, 20000, 21000);

        SvVarData var4 = createDup(tester.nextVarId(), chromosome, 30000, 31000);

        SvVarData var5 = createTestSv(tester.nextVarId(), chromosome, chromosome, 40000,45000, -1, -1, INV,  2.5);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        driverAnnotator.clearResults(false);
        driverAnnotator.annotateSVs(tester.SampleId, tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, driverAnnotator.getDriverGeneDataList().size());
        dgData = driverAnnotator.getDriverGeneDataList().get(0);

        assertEquals(2, dgData.getEvents().size());

        List<SvCluster> driverClusters = dgData.getEvents().stream().map(x -> x.getCluster()).collect(Collectors.toList());

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var2));
        assertTrue(cluster != null);
        assertTrue(driverClusters.contains(cluster));

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var5));
        assertTrue(cluster != null);
        assertTrue(driverClusters.contains(cluster));

        tester.clearClustersAndSVs();

        // on upper side: BND cancelled by opposing BND, INV partly cancelled
        var1 = createDel(tester.nextVarId(), chromosome, 500, 600);

        // BNDs are too far away to merge with the INV
        var2 = createTestSv(tester.nextVarId(), chromosome, "2", 6300000,1000, 1, 1, BND,  5);

        var3 = createTestSv(tester.nextVarId(), chromosome, "3", 6250000,1000, 1, 1, BND,  15);

        var4 = createTestSv(tester.nextVarId(), chromosome, chromosome, 200000,201000, -1, -1, INV,  5);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        driverAnnotator.clearResults(false);
        driverAnnotator.annotateSVs(tester.SampleId, tester.Analyser.getState().getChrBreakendMap());

        assertEquals(1, driverAnnotator.getDriverGeneDataList().size());
        dgData = driverAnnotator.getDriverGeneDataList().get(0);

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var3));

        assertEquals(1, dgData.getEvents().size());
        assertEquals(cluster, dgData.getEvents().get(0).getCluster());
    }

}
