package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createTestSv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import com.hartwig.hmftools.linx.utils.LinxTester;

// tests on chain functions and clusters not involving varied ploidy

public class ChainingSimpleTest
{
    @Test
    public void testFullyAssembledChain()
    {
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createDel(0, "1", 100,200);
        final SvVarData var2 = createDel(1, "1", 300,400);
        final SvVarData var3 = createDel(2, "1", 500,600);
        final SvVarData var4 = createDel(3, "1", 700,800);

        var1.setAssemblyData(false, "asmb12");
        var2.setAssemblyData(true, "asmb12");
        var2.setAssemblyData(false, "asmb23");
        var3.setAssemblyData(true, "asmb23");
        var3.setAssemblyData(false, "asmb34");
        var4.setAssemblyData(true, "asmb34");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(3, chain.getLinkCount());
        assertEquals(3, cluster.getAssemblyLinkedPairs().size());
    }

    @Test
    public void testChainNotClosed()
    {
        // 2 SVs which could link on both ends
        LinxTester tester = new LinxTester();
        // tester.logVerbose(true);

        final SvVarData var1 = createInv(1, "1", 100,200, -1);
        final SvVarData var2 = createDel(2, "1", 300,400);
        final SvVarData var3 = createDel(3, "1", 500,600);
        final SvVarData var4 = createInv(4, "1", 700,800, 1);

        // add them out of order which will require partial chain reconciliation
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(3, chain.getLinkCount());
    }

    @Test
    public void testPartiallyAssembledChain()
    {
        LinxTester tester = new LinxTester();
        // tester.logVerbose(true);

        final SvVarData var0 = createDel(0, "1", 100,200);
        final SvVarData var1 = createDel(1, "1", 300,400);
        final SvVarData var2 = createDel(2, "1", 500,600);
        final SvVarData var3 = createDel(3, "1", 700,800);

        var1.setAssemblyData(false, "asmb23");
        var2.setAssemblyData(true, "asmb23");

        // add them out of order which will require partial chain reconciliation
        tester.AllVariants.add(var0);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertEquals(1, cluster.getChains().size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(3, chain.getLinkCount());
        assertEquals(4, chain.getSvCount());
    }

    @Test
    public void testSkipSpanningAssemblies()
    {
        LinxTester tester = new LinxTester();

        final SvVarData var0 = createDel(0, "1", 100,200);
        final SvVarData var1 = createDel(1, "1", 300,400);
        final SvVarData var2 = createDel(2, "1", 500,600);
        final SvVarData var3 = createDel(3, "1", 700,800);
        final SvVarData var4 = createDel(4, "1", 900,1000);

        var0.setAssemblyData(false, "asmb01");
        var1.setAssemblyData(true, "asmb01");
        var1.setAssemblyData(false, "asmb12");
        var2.setAssemblyData(true, "asmb12");
        var2.setAssemblyData(false, "asmb23");
        var3.setAssemblyData(true, "asmb23");
        var3.setAssemblyData(false, "asmb34");
        var4.setAssemblyData(true, "asmb34");

        // and some spanning assemblies which will be ignored
        var0.setAssemblyData(false, "asmb02");
        var2.setAssemblyData(true, "asmb02");
        var2.setAssemblyData(false, "asmb24");
        var4.setAssemblyData(true, "asmb24");

        // add them out of order which will require partial chain reconciliation
        tester.AllVariants.add(var0);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(!cluster.requiresReplication());
        assertEquals(1, cluster.getChains().size());

        assertEquals(1, var0.getAssembledLinkedPairs(false).size());
        assertEquals(1, var2.getAssembledLinkedPairs(true).size());
        assertEquals(1, var2.getAssembledLinkedPairs(false).size());
        assertEquals(1, var4.getAssembledLinkedPairs(true).size());

        final SvChain chain = cluster.getChains().get(0);

        assertEquals(4, chain.getLinkCount());
        assertEquals(5, chain.getSvCount());
    }

    @Test
    public void testZeroBaseInferredChains()
    {
        LinxTester tester = new LinxTester();

        // Configurator.setRootLevel(Level.DEBUG);

        final SvVarData var0 = createTestSv(1, "1", "1", 1000,2000, 1, -1, DEL,  1);

        final SvVarData var1 = createTestSv(1, "1", "0", 21000,-1, 1, 0, INF,  1);
        final SvVarData var2 = createTestSv(2, "1", "0", 21000,-1, -1, 0, SGL,  1);

        final SvVarData var3 = createTestSv(3, "1", "2", 40000,100, 1, -1, BND, 4);
        final SvVarData var4 = createTestSv(4, "1", "0", 40000,-1, -1, 1, INF, 4);
        final SvVarData var5 = createTestSv(5, "1", "2", 44000,200, -1, 1, BND, 4);
        final SvVarData var6 = createTestSv(6, "1", "0", 45000,-1, 1, 1, SGL, 4);

        tester.AllVariants.add(var0);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(3, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(var1, var2));
        assertTrue(cluster != null);

        assertTrue(cluster.getChains().isEmpty());

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var3, var4, var5, var6));
        assertTrue(cluster != null);

        assertEquals(1, cluster.getChains().size());
        assertEquals(3, cluster.getChains().get(0).getLinkCount());

        final LinkedPair pair = var4.getLinkedPair(true);
        assertTrue(pair != null);
        assertTrue(pair.hasVariant(var3));
    }
}
