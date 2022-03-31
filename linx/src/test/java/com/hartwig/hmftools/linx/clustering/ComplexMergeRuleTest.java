package com.hartwig.hmftools.linx.clustering;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.CONSEC_BREAKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.FOLDBACKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.LOH_CHAIN;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.OVERLAP_FOLDBACKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.SATELLITE_SGL;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.TI_JCN_MATCH;
import static com.hartwig.hmftools.linx.types.SvVarData.ASSEMBLY_TYPE_EQV;
import static com.hartwig.hmftools.linx.types.SvVarData.SGL_CENTRO_SATELLITE;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSgl;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createTestSv;

import static org.junit.Assert.assertEquals;

import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class ComplexMergeRuleTest
{
    @Test
    public void testFoldbackMerge()
    {
        LinxTester tester = new LinxTester();

        // 2 clusters with foldbacks on the same arm are merged
        SvVarData inv1 = createInv(tester.nextVarId(), "1", 100, 200, -1);
        tester.AllVariants.add(inv1);

        SvVarData inv2 = createInv(tester.nextVarId(), "1", 20000, 20100, 1);
        tester.AllVariants.add(inv2);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.getClusters().size());
        assertTrue(inv1.hasClusterReason(FOLDBACKS));
        assertTrue(inv2.hasClusterReason(FOLDBACKS));
    }

    @Test
    public void testLohResolvingClusterMerge()
    {
        // merge clusters based on their SVs being required to stop a SV chaining through an LOH event
        LinxTester tester = new LinxTester();

        SvVarData var1 = createDel(1, "1", 10000, 60000);
        SvVarData var2 = createDup(2, "1", 20000, 50000);
        SvVarData var3 = createBnd(3, "1", 30000, 1, "2", 10000, 1);
        SvVarData var4 = createBnd(4, "1", 40000, -1, "2", 20000, 1);
        SvVarData var5 = createDel(5, "1", 70000, 120000);
        SvVarData var6 = createDup(6, "1", 80000, 110000);
        SvVarData var7 = createBnd(7, "1", 90000, 1, "2", 80000, 1);
        SvVarData var8 = createBnd(8, "1", 100000, -1, "2", 90000, -1);
        SvVarData var9 = createDup(9, "2", 30000, 60000);

        // does run unto the DUP in the LOH but isn't clustered since is simple
        SvVarData var10 = createDel(10, "2", 40000, 50000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);
        tester.AllVariants.add(var8);
        tester.AllVariants.add(var9);
        tester.AllVariants.add(var10);

        List<LohEvent> lohData = tester.CnDataLoader.getLohData();

        lohData.add(new LohEvent("1", 10000, 20000, "DEL", "DUP", 1, var1.id(), var2.id()));

        lohData.add(new LohEvent("1", 50000, 60000,"DUP", "DEL", 1, var2.id(), var1.id()));

        lohData.add(new LohEvent("1", 110000, 120000, "DUP", "DEL", 1, var6.id(), var5.id()));

        lohData.add(new LohEvent("2", 20000, 30000,"BND", "DUP", 1, var4.id(), var9.id()));

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(2, tester.Analyser.getClusters().size());
        assertTrue(var10.getCluster().getSvCount() == 1);

        assertTrue(var2.hasClusterReason(LOH_CHAIN));
        assertTrue(var2.getClusterReason().contains(String.valueOf(var4.id())));
        assertTrue(var4.hasClusterReason(LOH_CHAIN));
        assertTrue(var4.getClusterReason().contains(String.valueOf(var2.id())));
    }

    @Test
    public void testConsistentBreakendOverlapMerge()
    {
        LinxTester tester = new LinxTester();

        List<SvVarData> allVariants = Lists.newArrayList();

        // a cluster has 3 consecutive breakends which span other unresolved SV breakends, which are then merged in
        SvVarData consec1 = createBnd(tester.nextVarId(), "1", 100000, 1, "2", 100, -1);
        allVariants.add(consec1);

        SvVarData consec2 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 200, 1);
        allVariants.add(consec2);

        SvVarData consec3 = createInv(tester.nextVarId(), "1", 30000, 101000, 1);
        allVariants.add(consec3);

        // create some SV in resolved clusters which will be ignored - a DEL with external TI , a simple DEL and a low-qual
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 10000, 1, "3", 200, 1);
        allVariants.add(var1);

        SvVarData var2 = createBnd(tester.nextVarId(), "1", 10100, -1, "3", 100, -1);
        allVariants.add(var2);

        SvVarData var3 = createDel(tester.nextVarId(), "1", 60000, 60100);
        allVariants.add(var3);

        // eqv breakend will be ignored
        SvVarData var4 = createSgl(tester.nextVarId(), "1", 80000, -1);
        var4.setAssemblyData(true, ASSEMBLY_TYPE_EQV);
        allVariants.add(var4);

        // now some SV which will be overlapped by the consecutive breakends
        SvVarData overlap1 = createBnd(tester.nextVarId(), "1", 20000, -1, "4", 200, 1);
        allVariants.add(overlap1);

        SvVarData overlap2 = createSgl(tester.nextVarId(), "1", 40000, -1);
        allVariants.add(overlap2);

        tester.AllVariants.addAll(allVariants);
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(4, tester.getAllClusters().size());

        SvCluster mainCluster = tester.findClusterWithSVs(Lists.newArrayList(consec1, consec2, consec3, overlap1, overlap2));

        if(mainCluster == null)
            assertTrue(false);

        assertTrue(overlap1.hasClusterReason(CONSEC_BREAKS));
        assertTrue(overlap2.hasClusterReason(CONSEC_BREAKS));
    }

    @Test
    public void testFoldbacksStraddlingBreakendMerge()
    {
        LinxTester tester = new LinxTester();

        // a cluster with 2 foldbacks and an unclustered breakend in between

        // foldbacks clustered due to the foldback rule
        SvVarData var1 = createInv(tester.nextVarId(), "1", 1000, 2000, -1);
        SvVarData var2 = createInv(tester.nextVarId(), "1", 100000, 101000, -1);

        // simple del not clustered
        SvVarData var3 = createDel(tester.nextVarId(), "1", 20000, 21000);

        SvVarData var4 = createBnd(tester.nextVarId(), "1", 40000, -1, "3", 100, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.isFoldback());
        assertTrue(var2.isFoldback());

        assertEquals(2, tester.getClusters().size());

        SvCluster mainCluster = tester.findClusterWithSVs(Lists.newArrayList(var1, var2, var4));

        assertNotNull(mainCluster);

        assertTrue(var4.hasClusterReason(OVERLAP_FOLDBACKS));
        assertTrue(mainCluster.hasClusterReason(OVERLAP_FOLDBACKS));

        // don't merge if both foldbacks face away
        tester.clearClustersAndSVs();

        var1 = createInv(tester.nextVarId(), "1", 1000, 2000, 1);

        // var2 as before

        var3 = createBnd(tester.nextVarId(), "1", 40000, -1, "3", 100, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.isFoldback());
        assertTrue(var2.isFoldback());

        assertEquals(2, tester.getClusters().size());
    }

    @Test
    public void testChainedFoldbackMerge()
    {
        LinxTester tester = new LinxTester();

        // Configurator.setRootLevel(Level.DEBUG);

        // cluster foldbacks on the same arm even if they're chained without assembly

        // foldbacks clustered due to the foldback rule
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 1000, -1, "2", 100, -1);
        SvVarData var2 = createBnd(tester.nextVarId(), "1", 2000, -1, "2", 200, 1);

        // second cluster is sufficiently complex and distant to not trigger other rules
        SvVarData var3 = createBnd(tester.nextVarId(), "1", 6101000, 1, "4", 100, -1);
        SvVarData var4 = createBnd(tester.nextVarId(), "1", 6102000, 1, "4", 200, 1);

        SvVarData var5 = createDup(tester.nextVarId(), "1", 6105000, 6110000);
        SvVarData var6 = createDup(tester.nextVarId(), "1", 6108000, 6112000);
        SvVarData var7 = createDel(tester.nextVarId(), "1", 6111000, 6120000);
        SvVarData var8 = createSgl(tester.nextVarId(), "1", 6115000, -1);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2));
        tester.AllVariants.addAll(Lists.newArrayList(var3, var4, var5, var6, var7, var8));

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.isChainedFoldback());
        assertTrue(var2.isChainedFoldback());
        assertTrue(var3.isChainedFoldback());
        assertTrue(var4.isChainedFoldback());

        assertEquals(1, tester.getClusters().size());
        assertTrue(var1.hasClusterReason(FOLDBACKS));
    }

    @Test
    public void testDistancePloidyLinkMatchMerge()
    {
        LinxTester tester = new LinxTester();

        // set a high major allele ploidy on the other chromatid to stop the MAP merge rule kicking in
        tester.setNonClusterAllelePloidies(4, 1);

        // 2 distance clusters which can form a distant ploidy-matching merge

        // foldbacks clustered due to the foldback rule
        SvVarData var0 = createTestSv(tester.nextVarId(), "1", "0", 500, -1, -1, -1, SGL, 2);
        SvVarData var1 = createTestSv(tester.nextVarId(), "1", "2", 1000, 100, 1, -1, BND, 2);
        SvVarData var2 = createTestSv(tester.nextVarId(), "1", "2", 2000, 2000, -1, 1, BND, 2);

        // simple del not clustered
        SvVarData var3 = createDel(tester.nextVarId(), "1", 20000, 21000);

        // pair of BNDs resolved to a simple type
        SvVarData var4 = createTestSv(tester.nextVarId(), "1", "3", 50000, 100, 1, -1, BND, 1);
        SvVarData var5 = createTestSv(tester.nextVarId(), "1", "3", 51000, 200, -1, 1, BND, 1);

        SvVarData var6 = createTestSv(tester.nextVarId(), "1", "1", 100000, 110000, 1, 1, INV, 2);

        tester.AllVariants.add(var0);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(3, tester.getClusters().size());

        assertTrue(tester.hasClusterWithSVs(Lists.newArrayList(var4, var5)));

        SvCluster mainCluster = tester.findClusterWithSVs(Lists.newArrayList(var0, var1, var2, var6));

        assertNotNull(mainCluster);

        assertTrue(var6.hasClusterReason(TI_JCN_MATCH));
        assertTrue(mainCluster.hasClusterReason(TI_JCN_MATCH));
    }

    @Test
    public void testSatelliteSglMerge()
    {
        LinxTester tester = new LinxTester();

        // foldbacks clustered due to the foldback rule
        SvVarData var0 = createTestSv(tester.nextVarId(), "1", "0", 500, -1, 1, 0, SGL,
                2, 0, 1, 0, 1, "", PASS, SGL_CENTRO_SATELLITE, "SAR");

        SvVarData var1 = createTestSv(tester.nextVarId(), "1", "0", 50000, -1, -1, 0, SGL,
                2, 0, 1, 0, 1, "", PASS, "", "SAR");

        tester.AllVariants.add(var0);
        tester.AllVariants.add(var1);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.getClusters().size());

        assertTrue(var0.hasClusterReason(SATELLITE_SGL));
        assertTrue(var1.hasClusterReason(SATELLITE_SGL));
    }
}
