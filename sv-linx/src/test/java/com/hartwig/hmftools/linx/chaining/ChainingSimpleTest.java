package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.chaining.SvChain.reconcileChains;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.types.SvCluster;
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
        // tester.logVerbose(true);

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
}
