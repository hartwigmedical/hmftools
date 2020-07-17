package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInf;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class LinkingTest
{
    @Test
    public void testLinkedPairs()
    {
        // test linked pair switching
        final SvVarData var1 = createDup(1, "1", 100, 200);
        final SvVarData var2 = createDup(2, "1", 300, 400);
        LinkedPair lp1 = LinkedPair.from(var1, var2, false, true);

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
        final SvVarData var3 = createDel(3, "1", 100, 200);
        final SvVarData var4 = createDel(4, "1", 210, 400);
        LinkedPair lp2 = LinkedPair.from(var3, var4, false, true);

        final SvVarData var5 = createDel(4, "1", 250, 400);

        LinkedPair lp3 = LinkedPair.from(var3, var5, false, true);

        lp3.sameVariants(lp2);
        lp3.hasLinkClash(lp2);
    }

    @Test
    public void testDeletionBridges()
    {
        LinxTester tester = new LinxTester();

        SvVarData var1 = createInv(tester.nextVarId(), "1", 100, 1000, 1);
        SvVarData var2 = createInf(tester.nextVarId(), "1", 1001, -1);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2));

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.getDBLink(false) != null);
        assertEquals(var1.getDBLink(false), var2.getDBLink(true));
        assertEquals(0, var1.getDBLink(false).length());

        tester.clearClustersAndSVs();

        var1 = createInv(tester.nextVarId(), "1", 100, 1000, 1);
        var2 = createInf(tester.nextVarId(), "1", 1000, -1);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2));

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.getDBLink(false) != null);
        assertEquals(var1.getDBLink(false), var2.getDBLink(true));
        assertEquals(-1, var1.getDBLink(false).length());


    }

}
