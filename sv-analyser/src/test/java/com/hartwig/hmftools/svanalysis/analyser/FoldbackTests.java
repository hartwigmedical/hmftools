package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createSgl;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_INT_TI;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Test;

public class FoldbackTests
{

    @Test
    public void testInvFoldbacks()
    {
        SvTestHelper tester = new SvTestHelper();

        // create foldbacks from a simple inversion

        SvVarData var1 = createInv(tester.nextVarId(), "1", 100, 200, 1);

        SvCluster cluster = new SvCluster(0);
        cluster.addVariant(var1);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();
        tester.Analyser.markFoldbacks();

        assertEquals(var1.getFoldbackLink(true), var1.id());
        assertEquals(var1.getFoldbackLink(false), var1.id());
        assertEquals(var1.getFoldbackLen(true), 100);

        // now test with a deletion bridge obscuring the backmost breakend

        tester.clearClustersAndSVs();

        SvVarData dup = createDup(tester.nextVarId(), "1", 190, 400);
        cluster.addVariant(dup);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();
        tester.Analyser.markFoldbacks();

        assertEquals(var1.getFoldbackLink(true), var1.id());
        assertEquals(var1.getFoldbackLink(false), var1.id());
        assertEquals(var1.getFoldbackLen(true), 100);

        // now invalidate the foldback by putting a deletion bridge at the start
        tester.clearClustersAndSVs();

        SvVarData var2 = createInv(tester.nextVarId(), "1", 100, 200, -1);

        SvCluster cluster2 = new SvCluster(1);
        cluster2.addVariant(var2);

        SvVarData del = createDel(tester.nextVarId(), "1", 190, 400);
        cluster2.addVariant(del);

        tester.addClusterAndSVs(cluster2);
        tester.preClusteringInit();

        SvLinkedPair dbPair = var2.getDBLink(false);
        assertTrue(dbPair != null);
        assertTrue(dbPair.hasBreakend(var2, false));
        assertTrue(dbPair.hasBreakend(del, true));

        tester.Analyser.markFoldbacks();

        assertEquals(var2.getFoldbackLink(true), "");
        assertEquals(var2.getFoldbackLink(false), "");

        // finally an INV with DBs at both ends, making it invalid
        tester.clearClustersAndSVs();

        SvVarData var3 = createInv(tester.nextVarId(), "1", 100, 200, 1);

        SvCluster cluster3 = new SvCluster(2);
        cluster3.addVariant(var3);

        SvVarData sgl = createSgl(tester.nextVarId(), "1", 90, -1, false);
        cluster3.addVariant(sgl);

        SvVarData sgl2 = createSgl(tester.nextVarId(), "1", 190, -1, false);
        cluster3.addVariant(sgl2);

        tester.addClusterAndSVs(cluster3);
        tester.preClusteringInit();

        dbPair = var3.getDBLink(true);
        assertTrue(dbPair != null);
        assertTrue(dbPair.hasBreakend(var3, true));
        assertTrue(dbPair.hasBreakend(sgl, true));

        dbPair = var3.getDBLink(false);
        assertTrue(dbPair != null);
        assertTrue(dbPair.hasBreakend(var3, false));
        assertTrue(dbPair.hasBreakend(sgl2, true));

        tester.Analyser.markFoldbacks();

        assertEquals(var3.getFoldbackLink(true), "");
        assertEquals(var3.getFoldbackLink(false), "");
    }

}