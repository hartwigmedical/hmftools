package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.linx.chaining.SvChain.CHAIN_ASSEMBLY_LINK_COUNT;
import static com.hartwig.hmftools.linx.chaining.SvChain.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.linx.chaining.SvChain.checkIsValid;
import static com.hartwig.hmftools.linx.chaining.SvChain.reconcileChains;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

// tests on chain functions and clusters not involving varied ploidy

public class ChainingSimpleTest
{
    @Test
    public void testChainRoutines()
    {
        // create a chain out of simple DELs and test the various chaining features
        final SvVarData var1 = createDel("1", "1", 1100, 1200);
        final SvVarData var2 = createDel("2", "1", 1300, 1400);
        SvLinkedPair lp1 = new SvLinkedPair(var1, var2, LINK_TYPE_TI, false, true);
        lp1.setIsAssembled();

        // test adding linked pairs of various orientations to the start and end of a chain
        SvChain chain = new SvChain(0);

        chain.addLink(lp1, true);

        assertTrue(chain.firstLinkOpenOnStart());
        assertTrue(!chain.lastLinkOpenOnStart());

        final SvVarData var3 = createDel("3", "1", 1500, 1600);
        SvLinkedPair lp2 = new SvLinkedPair(var3, var2, LINK_TYPE_TI, true, false);
        lp2.setIsAssembled();

        assertFalse(chain.canAddLinkedPairToStart(lp2));
        assertTrue(chain.canAddLinkedPairToEnd(lp2));

        chain.addLink(lp2, false);

        assertEquals(chain.getFirstSV(), var1);
        assertEquals(chain.getLastSV(), var3);
        assertEquals(chain.getSvList().get(1), var2);

        final SvVarData var4 = createDel("4", "1", 900, 1000);
        SvLinkedPair lp3 = new SvLinkedPair(var1, var4, LINK_TYPE_TI, true, false);

        assertTrue(chain.canAddLinkedPairToStart(lp3));
        assertFalse(chain.canAddLinkedPairToEnd(lp3));
        chain.addLink(lp3, true);

        assertTrue(checkIsValid(chain));

        assertEquals(chain.getFirstSV(), var4);
        assertEquals(chain.getLastSV(), var3);

        // test a potentially closing link
        final SvVarData var5 = createDup("5", "1", 800, 1700);
        SvLinkedPair lp4 = new SvLinkedPair(var5, var4, LINK_TYPE_TI, true, true);

        assertTrue(chain.canAddLinkedPairToStart(lp4));
        chain.addLink(lp4, true);

        assertEquals(chain.getAssemblyLinkCount(), 2);

        assertTrue(checkIsValid(chain));

        SvLinkedPair lp5 = new SvLinkedPair(var5, var3, LINK_TYPE_TI, false, false);

        assertTrue(chain.canAddLinkedPairToEnd(lp5));
        assertTrue(chain.linkWouldCloseChain(lp5));

        chain.closeChain("CLOSE", 0);
        assertTrue(chain.isClosedLoop());
        assertEquals(5, chain.getLinkCount());

        // tests paths through the chain from various points
        int[] chainData = chain.breakendsAreChained(var4, false, var3, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 3);
        assertEquals(chainData[CHAIN_ASSEMBLY_LINK_COUNT], 2);

        // check works in the other direction
        chainData = chain.breakendsAreChained(var3, true, var4, false);
        assertEquals(chainData[CHAIN_LINK_COUNT], 3);
        assertEquals(chainData[CHAIN_ASSEMBLY_LINK_COUNT], 2);

        // check a single link
        chainData = chain.breakendsAreChained(var1, false, var2, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 1);
        assertEquals(chainData[CHAIN_ASSEMBLY_LINK_COUNT], 1);

        // check breakends facing the wrong way
        chainData = chain.breakendsAreChained(var1, false, var2, false);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);

        chainData = chain.breakendsAreChained(var1, true, var2, false);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);

        chainData = chain.breakendsAreChained(var1, true, var2, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);

        // check no link
        chainData = chain.breakendsAreChained(var5, false, var1, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);
    }

    @Test
    public void testChainReconciliation()
    {
        SvVarData var1 = createDel("1", "1", 100, 150);
        SvVarData var2 = createDel("2", "1", 200, 250);
        SvVarData var3 = createDel("3", "1", 300, 350);
        SvVarData var4 = createDel("4", "1", 400, 450);
        SvVarData var5 = createDel("5", "1", 500, 550);
        SvVarData var6 = createDel("6", "1", 600, 650);
        SvVarData var7 = createDel("7", "1", 700, 750);

        SvChain chain1 = new SvChain(1);
        chain1.addLink(SvLinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain1.addLink(SvLinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        SvChain chain2 = new SvChain(2);
        chain2.addLink(SvLinkedPair.from(var3.getBreakend(false), var4.getBreakend(true)), true);
        chain2.addLink(SvLinkedPair.from(var4.getBreakend(false), var5.getBreakend(true)), false);

        List<SvChain> chainList = Lists.newArrayList(chain1, chain2);

        reconcileChains(chainList);

        assertTrue(chainList.size() == 1);
        assertEquals(4, chain1.getLinkCount());
        assertEquals(5, chain1.getSvCount());
        assertTrue(checkIsValid(chain1));

        // reset and try again with 3 chains
        chain1 = new SvChain(1);
        chain1.addLink(SvLinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain1.addLink(SvLinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        SvChain chain3 = new SvChain(3);
        chain3.addLink(SvLinkedPair.from(var5.getBreakend(false), var6.getBreakend(true)), true);
        chain3.addLink(SvLinkedPair.from(var6.getBreakend(false), var7.getBreakend(true)), false);

        chainList = Lists.newArrayList(chain3, chain1, chain2);

        reconcileChains(chainList);

        assertTrue(chainList.size() == 1);
        assertEquals(6, chain3.getLinkCount());
        assertEquals(7, chain3.getSvCount());
        assertTrue(checkIsValid(chain3));

        // check cannot merge chains on the inner starting SVs
        chain1 = new SvChain(1);
        chain1.addLink(SvLinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain1.addLink(SvLinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        chain2 = new SvChain(2);
        chain2.addLink(SvLinkedPair.from(var2.getBreakend(false), var4.getBreakend(true)), true);
        chain2.addLink(SvLinkedPair.from(var4.getBreakend(false), var5.getBreakend(true)), false);

        chainList = Lists.newArrayList(chain1, chain2);
        reconcileChains(chainList);

        assertTrue(chainList.size() == 2);

        chainList = Lists.newArrayList(chain2, chain1);
        reconcileChains(chainList);

        assertTrue(chainList.size() == 2);

    }


    @Test
    public void testChainSplittingByFoldback()
    {
        SvTestHelper tester = new SvTestHelper();

        final SvVarData var1 = createDel("1", "1", 300, 400);
        final SvVarData var2 = createDel("2", "1", 500, 600);
        final SvVarData var3 = createDel("3", "1", 700, 800);

        // test adding linked pairs of various orientations to the start and end of a chain
        SvChain chain = new SvChain(0);

        chain.addLink(SvLinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain.addLink(SvLinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        assertTrue(checkIsValid(chain));

        // now add a foldback which will split and replicate the chain
        SvVarData var4 = createInv("4", "1", 900, 1000, 1);

        chain.foldbackChainOnLink(
                SvLinkedPair.from(var4.getBreakend(true), var3.getBreakend(false)),
                SvLinkedPair.from(var4.getBreakend(false), var3.getBreakend(false)));

        assertEquals(6, chain.getLinkCount());
        assertEquals(4, chain.getSvCount());
        assertTrue(checkIsValid(chain));

        // again but with a foldback connected at the start
        chain = new SvChain(0);

        chain.addLink(SvLinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain.addLink(SvLinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        assertTrue(checkIsValid(chain));

        // now add a foldback which will split and replicate the chain
        var4 = createInv("4", "1", 100, 200, -1);

        chain.foldbackChainOnLink(
                SvLinkedPair.from(var1.getBreakend(true), var4.getBreakend(true)),
                SvLinkedPair.from(var4.getBreakend(false), var1.getBreakend(true)));

        assertEquals(6, chain.getLinkCount());
        assertEquals(4, chain.getSvCount());
        assertTrue(checkIsValid(chain));

    }

    @Test
    public void testChainSplittingByComplexDup()
    {
        SvTestHelper tester = new SvTestHelper();

        final SvVarData var1 = createDel("1", "1", 300, 400);
        final SvVarData var2 = createDel("2", "1", 500, 600);
        final SvVarData var3 = createDel("3", "1", 700, 800);

        // test adding linked pairs of various orientations to the start and end of a chain
        SvChain chain = new SvChain(0);

        chain.addLink(SvLinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain.addLink(SvLinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        assertTrue(checkIsValid(chain));

        // now add a foldback which will split and replicate the chain
        SvVarData var4 = createDup("4", "1", 200, 900);

        chain.duplicateChainOnLink(
                SvLinkedPair.from(var4.getBreakend(true), var1.getBreakend(true)),
                SvLinkedPair.from(var4.getBreakend(false), var3.getBreakend(false)));

        assertEquals(6, chain.getLinkCount());
        assertEquals(4, chain.getSvCount());
        assertTrue(checkIsValid(chain));
    }

    @Test
    public void testFullyAssembledChain()
    {
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        final SvVarData var1 = createDel("0", "1", 100,200);
        final SvVarData var2 = createDel("1", "1", 300,400);
        final SvVarData var3 = createDel("2", "1", 500,600);
        final SvVarData var4 = createDel("3", "1", 700,800);

        var1.setAssemblyData(false, "asmb12");
        var2.setAssemblyData(true, "asmb12");
        var2.setAssemblyData(false, "asmb23");
        var3.setAssemblyData(true, "asmb23");
        var3.setAssemblyData(false, "asmb34");
        var4.setAssemblyData(true, "asmb34");

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
        assertEquals(3, cluster.getAssemblyLinkedPairs().size());
    }

    @Test
    public void testChainNotClosed()
    {
        // 2 SVs which could link on both ends
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        final SvVarData var1 = createInv("0", "1", 100,200, -1);
        final SvVarData var2 = createDel("1", "1", 300,400);
        final SvVarData var3 = createDel("2", "1", 500,600);
        final SvVarData var4 = createInv("3", "1", 700,800, 1);

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
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        final SvVarData var0 = createDel("0", "1", 100,200);
        final SvVarData var1 = createDel("1", "1", 300,400);
        final SvVarData var2 = createDel("2", "1", 500,600);
        final SvVarData var3 = createDel("3", "1", 700,800);

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
    }
}
