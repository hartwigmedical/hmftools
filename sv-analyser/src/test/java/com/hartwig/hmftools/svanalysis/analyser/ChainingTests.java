package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_ASSEMBLY_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_DB;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Test;

public class ChainingTests
{
    @Test
    public void testLinkedPairs()
    {
        // test linked pair switching
        final SvVarData var1 = createDup("1", "1", 100, 200);
        final SvVarData var2 = createDup("2", "1", 300, 400);
        SvLinkedPair lp1 = new SvLinkedPair(var1, var2, LINK_TYPE_TI, false, true);

        assertEquals(lp1.first(), var1);
        assertEquals(lp1.second(), var2);
        assertEquals(lp1.firstLinkOnStart(), false);
        assertEquals(lp1.secondLinkOnStart(), true);
        assertEquals(lp1.firstUnlinkedOnStart(), true);
        assertEquals(lp1.secondUnlinkedOnStart(), false);

        lp1.switchSVs();
        assertEquals(lp1.first(), var2);
        assertEquals(lp1.second(), var1);
        assertEquals(lp1.firstLinkOnStart(), true);
        assertEquals(lp1.secondLinkOnStart(), false);

        // test short TIs converted to DBs
        final SvVarData var3 = createDel("3", "1", 100, 200);
        final SvVarData var4 = createDel("4", "1", 210, 400);
        SvLinkedPair lp2 = new SvLinkedPair(var3, var4, LINK_TYPE_TI, false, true);
        assertEquals(lp2.linkType(), LINK_TYPE_DB);
        assertEquals(lp2.length(), -10);

        final SvVarData var5 = createDel("4", "1", 250, 400);

        SvLinkedPair lp3 = new SvLinkedPair(var3, var5, LINK_TYPE_TI, false, true);

        lp3.sameVariants(lp2);
        lp3.hasLinkClash(lp2);
    }

    @Test
    public void testBasicChainBuilding()
    {
        // create a chain out of simple DELs and test the various chaining features
        final SvVarData var1 = createDel("1", "1", 1100, 1200);
        final SvVarData var2 = createDel("2", "1", 1300, 1400);
        SvLinkedPair lp1 = new SvLinkedPair(var1, var2, LINK_TYPE_TI, false, true);
        lp1.setIsInferred(false);

        // test adding linked pairs of various orientations to the start and end of a chain
        SvChain chain = new SvChain(0);

        chain.addLink(lp1, true);

        assertTrue(chain.firstLinkOpenOnStart());
        assertTrue(!chain.lastLinkOpenOnStart());

        final SvVarData var3 = createDel("3", "1", 1500, 1600);
        SvLinkedPair lp2 = new SvLinkedPair(var3, var2, LINK_TYPE_TI, true, false);
        lp2.setIsInferred(false);

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

        assertEquals(chain.getFirstSV(), var4);
        assertEquals(chain.getLastSV(), var3);
        assertEquals(chain.getSvList().get(1), var1);
        assertEquals(chain.getSvList().get(2), var2);

        // test a potentially closing link
        final SvVarData var5 = createDup("5", "1", 800, 1700);
        SvLinkedPair lp4 = new SvLinkedPair(var5, var4, LINK_TYPE_TI, true, true);

        assertTrue(chain.canAddLinkedPairToStart(lp4));
        chain.addLink(lp4, true);

        assertEquals(chain.getAssemblyLinkCount(), 2);

        SvLinkedPair lp5 = new SvLinkedPair(var5, var3, LINK_TYPE_TI, false, false);

        assertTrue(chain.canAddLinkedPairToEnd(lp5));
        assertTrue(chain.linkWouldCloseChain(lp5));

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


}
