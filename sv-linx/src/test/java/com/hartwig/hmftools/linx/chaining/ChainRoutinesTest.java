package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.linx.chaining.ChainUtils.CHAIN_ASSEMBLY_LINK_COUNT;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.breakendsAreChained;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.duplicateChainOnLink;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.foldbackChainOnChain;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.foldbackChainOnLink;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.reconcileChains;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.reverseSectionOnBreakend;
import static com.hartwig.hmftools.linx.chaining.SvChain.checkIsValid;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

public class ChainRoutinesTest
{
    @Test
    public void testChainRoutines()
    {
        // create a chain out of simple DELs and test the various chaining features
        final SvVarData var1 = createDel(1, "1", 1100, 1200);
        final SvVarData var2 = createDel(2, "1", 1300, 1400);
        LinkedPair lp1 = LinkedPair.from(var1, var2, false, true);
        lp1.setIsAssembled();

        // test adding linked pairs of various orientations to the start and end of a chain
        SvChain chain = new SvChain(0);

        chain.addLink(lp1, true);

        assertTrue(chain.firstLinkOpenOnStart());
        assertTrue(!chain.lastLinkOpenOnStart());

        final SvVarData var3 = createDel(3, "1", 1500, 1600);
        LinkedPair lp2 = LinkedPair.from(var3, var2, true, false);
        lp2.setIsAssembled();

        assertFalse(chain.canAddLinkedPairToStart(lp2));
        assertTrue(chain.canAddLinkedPairToEnd(lp2));

        chain.addLink(lp2, false);

        assertEquals(chain.getFirstSV(), var1);
        assertEquals(chain.getLastSV(), var3);
        assertEquals(chain.getSvList().get(1), var2);

        final SvVarData var4 = createDel(4, "1", 900, 1000);
        LinkedPair lp3 = LinkedPair.from(var1, var4, true, false);

        assertTrue(chain.canAddLinkedPairToStart(lp3));
        assertFalse(chain.canAddLinkedPairToEnd(lp3));
        chain.addLink(lp3, true);

        assertTrue(checkIsValid(chain));

        assertEquals(chain.getFirstSV(), var4);
        assertEquals(chain.getLastSV(), var3);

        // test a potentially closing link
        final SvVarData var5 = createDup(5, "1", 800, 1700);
        LinkedPair lp4 = LinkedPair.from(var5, var4, true, true);

        assertTrue(chain.canAddLinkedPairToStart(lp4));
        chain.addLink(lp4, true);

        assertEquals(chain.getAssemblyLinkCount(), 2);

        assertTrue(checkIsValid(chain));

        LinkedPair lp5 = LinkedPair.from(var5, var3, false, false);

        assertTrue(chain.canAddLinkedPairToEnd(lp5));
        assertTrue(chain.linkWouldCloseChain(lp5));

        chain.closeChain("CLOSE", 0);
        assertTrue(chain.isClosedLoop());
        assertEquals(5, chain.getLinkCount());

        // tests paths through the chain from various points
        int[] chainData = breakendsAreChained(chain, var4, false, var3, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 3);
        assertEquals(chainData[CHAIN_ASSEMBLY_LINK_COUNT], 2);

        // check works in the other direction
        chainData = breakendsAreChained(chain, var3, true, var4, false);
        assertEquals(chainData[CHAIN_LINK_COUNT], 3);
        assertEquals(chainData[CHAIN_ASSEMBLY_LINK_COUNT], 2);

        // check a single link
        chainData = breakendsAreChained(chain, var1, false, var2, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 1);
        assertEquals(chainData[CHAIN_ASSEMBLY_LINK_COUNT], 1);

        // check breakends facing the wrong way
        chainData = breakendsAreChained(chain, var1, false, var2, false);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);

        chainData = breakendsAreChained(chain, var1, true, var2, false);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);

        chainData = breakendsAreChained(chain, var1, true, var2, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);

        // check no link
        chainData = breakendsAreChained(chain, var5, false, var1, true);
        assertEquals(chainData[CHAIN_LINK_COUNT], 0);
    }

    @Test
    public void testChainReconciliation()
    {
        SvVarData var1 = createDel(1, "1", 100, 150);
        SvVarData var2 = createDel(2, "1", 200, 250);
        SvVarData var3 = createDel(3, "1", 300, 350);
        SvVarData var4 = createDel(4, "1", 400, 450);
        SvVarData var5 = createDel(5, "1", 500, 550);
        SvVarData var6 = createDel(6, "1", 600, 650);
        SvVarData var7 = createDel(7, "1", 700, 750);

        SvChain chain1 = new SvChain(1);
        chain1.addLink(LinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain1.addLink(LinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        SvChain chain2 = new SvChain(2);
        chain2.addLink(LinkedPair.from(var3.getBreakend(false), var4.getBreakend(true)), true);
        chain2.addLink(LinkedPair.from(var4.getBreakend(false), var5.getBreakend(true)), false);

        List<SvChain> chainList = Lists.newArrayList(chain1, chain2);

        reconcileChains(chainList);

        assertTrue(chainList.size() == 1);
        assertEquals(4, chain1.getLinkCount());
        assertEquals(5, chain1.getSvCount());
        assertTrue(checkIsValid(chain1));

        // reset and try again with 3 chains
        chain1 = new SvChain(1);
        chain1.addLink(LinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain1.addLink(LinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        SvChain chain3 = new SvChain(3);
        chain3.addLink(LinkedPair.from(var5.getBreakend(false), var6.getBreakend(true)), true);
        chain3.addLink(LinkedPair.from(var6.getBreakend(false), var7.getBreakend(true)), false);

        chainList = Lists.newArrayList(chain3, chain1, chain2);

        reconcileChains(chainList);

        assertTrue(chainList.size() == 1);
        assertEquals(6, chain3.getLinkCount());
        assertEquals(7, chain3.getSvCount());
        assertTrue(checkIsValid(chain3));

        // check cannot merge chains on the inner starting SVs
        chain1 = new SvChain(1);
        chain1.addLink(LinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain1.addLink(LinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        chain2 = new SvChain(2);
        chain2.addLink(LinkedPair.from(var2.getBreakend(false), var4.getBreakend(true)), true);
        chain2.addLink(LinkedPair.from(var4.getBreakend(false), var5.getBreakend(true)), false);

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
        final SvVarData var1 = createDel(1, "1", 300, 400);
        final SvVarData var2 = createDel(2, "1", 500, 600);
        final SvVarData var3 = createDel(3, "1", 700, 800);

        SvChain chain = new SvChain(0);

        chain.addLink(LinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain.addLink(LinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        assertTrue(checkIsValid(chain));

        // now add a foldback which will split and replicate the chain
        SvVarData var4 = createInv(4, "1", 900, 1000, 1);

        foldbackChainOnLink(chain,
                LinkedPair.from(var4.getBreakend(true), var3.getBreakend(false)),
                LinkedPair.from(var4.getBreakend(false), var3.getBreakend(false)));

        assertEquals(6, chain.getLinkCount());
        assertEquals(4, chain.getSvCount());
        assertTrue(checkIsValid(chain));
        assertEquals(chain.getOpenBreakend(true), chain.getOpenBreakend(false));

        // again but with a foldback connected at the start
        chain = new SvChain(0);

        chain.addLink(LinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain.addLink(LinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        assertTrue(checkIsValid(chain));

        // now add a foldback which will split and replicate the chain
        var4 = createInv(4, "1", 100, 200, -1);

        foldbackChainOnLink(chain,
                LinkedPair.from(var1.getBreakend(true), var4.getBreakend(true)),
                LinkedPair.from(var4.getBreakend(false), var1.getBreakend(true)));

        assertEquals(6, chain.getLinkCount());
        assertEquals(4, chain.getSvCount());
        assertTrue(checkIsValid(chain));
        assertEquals(chain.getOpenBreakend(true), chain.getOpenBreakend(false));
    }

    @Test
    public void testChainSplittingByFoldbackChain()
    {
        final SvVarData var1 = createDel(1, "1", 300, 400);
        final SvVarData var2 = createDel(2, "1", 500, 600);
        final SvVarData var3 = createDel(3, "1", 700, 800);

        SvChain chain = new SvChain(0);

        chain.addLink(LinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain.addLink(LinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        assertTrue(checkIsValid(chain));

        // now add a foldback chain which will split and replicate the chain
        SvChain fbChain = new SvChain(1);

        SvVarData var4 = createBnd(4, "1", 900, 1, "2", 100, -1);
        SvVarData var5 = createDel(5, "2", 500, 600);
        SvVarData var6 = createBnd(6, "1", 1000, -1, "2", 1000, 1);

        fbChain.addLink(LinkedPair.from(var4.getBreakend(false), var5.getBreakend(true)), true);
        fbChain.addLink(LinkedPair.from(var5.getBreakend(false), var6.getBreakend(false)), false);

        foldbackChainOnChain(chain,
                fbChain,
                LinkedPair.from(var4.getBreakend(true), var3.getBreakend(false)),
                LinkedPair.from(var6.getBreakend(true), var3.getBreakend(false)));

        assertEquals(8, chain.getLinkCount());
        assertEquals(6, chain.getSvCount());
        assertTrue(checkIsValid(chain));
        assertEquals(chain.getOpenBreakend(true), chain.getOpenBreakend(false));

        // again but with a foldback connected at the start
        chain = new SvChain(0);

        chain.addLink(LinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain.addLink(LinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        assertTrue(checkIsValid(chain));

        // now add a foldback which will split and replicate the chain
        fbChain = new SvChain(1);

        var4 = createBnd(4, "1", 100, 1, "2", 100, -1);
        var5 = createDel(5, "2", 500, 600);
        var6 = createBnd(6, "1", 200, -1, "2", 1000, 1);

        fbChain.addLink(LinkedPair.from(var4.getBreakend(false), var5.getBreakend(true)), true);
        fbChain.addLink(LinkedPair.from(var5.getBreakend(false), var6.getBreakend(false)), false);

        foldbackChainOnChain(chain,
                fbChain,
                LinkedPair.from(var1.getBreakend(true), var4.getBreakend(true)),
                LinkedPair.from(var1.getBreakend(true), var6.getBreakend(true)));

        assertEquals(8, chain.getLinkCount());
        assertEquals(6, chain.getSvCount());
        assertTrue(checkIsValid(chain));

        assertEquals(chain.getOpenBreakend(true), chain.getOpenBreakend(false));
    }

    @Test
    public void testChainSplittingByComplexDup()
    {
        final SvVarData var1 = createDel(1, "1", 300, 400);
        final SvVarData var2 = createDel(2, "1", 500, 600);
        final SvVarData var3 = createDel(3, "1", 700, 800);

        SvChain chain = new SvChain(0);

        chain.addLink(LinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain.addLink(LinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);

        assertTrue(checkIsValid(chain));

        // now add a complex dup which will split and replicate the chain
        SvVarData var4 = createDup(4, "1", 200, 900);

        duplicateChainOnLink(chain,
                LinkedPair.from(var4.getBreakend(true), var1.getBreakend(true)),
                LinkedPair.from(var4.getBreakend(false), var3.getBreakend(false)));

        assertEquals(6, chain.getLinkCount());
        assertEquals(4, chain.getSvCount());
        assertTrue(checkIsValid(chain));
    }

    @Test
    public void testChainSectionReversal()
    {
        SvVarData var1 = createInv(1, "1", 100, 200, -1);
        SvVarData var2 = createInv(2, "1", 300, 400, 1);

        SvChain chain = new SvChain(0);

        // test with a single-link chain
        chain.addLink(LinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);

        reverseSectionOnBreakend(chain, var1.getBreakend(false));
        assertEquals(var1.getBreakend(false), chain.getOpenBreakend(true));
        assertEquals(var2.getBreakend(false), chain.getOpenBreakend(false));

        reverseSectionOnBreakend(chain, var2.getBreakend(true));
        assertEquals(var1.getBreakend(false), chain.getOpenBreakend(true));
        assertEquals(var2.getBreakend(true), chain.getOpenBreakend(false));

        // now with a longer chain reversed in the middle
        var1 = createDel(1, "1", 100, 200);
        var2 = createDel(2, "1", 300, 400);
        final SvVarData var3 = createDel(3, "1", 500, 600);
        final SvVarData var4 = createDel(4, "1", 700, 800);
        final SvVarData var5 = createDel(5, "1", 900, 1000);

        chain = new SvChain(0);

        chain.addLink(LinkedPair.from(var1.getBreakend(false), var2.getBreakend(true)), true);
        chain.addLink(LinkedPair.from(var2.getBreakend(false), var3.getBreakend(true)), false);
        chain.addLink(LinkedPair.from(var3.getBreakend(false), var4.getBreakend(true)), false);
        chain.addLink(LinkedPair.from(var4.getBreakend(false), var5.getBreakend(true)), false);


        reverseSectionOnBreakend(chain, var3.getBreakend(false));
        assertTrue(checkIsValid(chain));
        assertEquals(var3.getBreakend(false), chain.getOpenBreakend(true));

        reverseSectionOnBreakend(chain, var1.getBreakend(true));
        assertTrue(checkIsValid(chain));
        assertEquals(var1.getBreakend(true), chain.getOpenBreakend(true));

        reverseSectionOnBreakend(chain, var2.getBreakend(true));
        assertTrue(checkIsValid(chain));
        assertEquals(var2.getBreakend(true), chain.getOpenBreakend(false));

        reverseSectionOnBreakend(chain, var5.getBreakend(false));
        assertTrue(checkIsValid(chain));
        assertEquals(var5.getBreakend(false), chain.getOpenBreakend(false));

    }
}
