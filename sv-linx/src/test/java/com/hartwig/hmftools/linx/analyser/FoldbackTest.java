package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.linx.utils.SvTestRoutines.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestRoutines.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestRoutines.createInv;
import static com.hartwig.hmftools.linx.utils.SvTestRoutines.createSgl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.analysis.FoldbackFinder;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import com.hartwig.hmftools.linx.utils.LinxTester;

public class FoldbackTest
{

    @Test
    public void testInvFoldbacks()
    {
        LinxTester tester = new LinxTester();

        // create foldbacks from a simple inversion

        SvVarData var1 = createInv(tester.nextVarId(), "1", 100, 200, 1);

        SvCluster cluster1 = new SvCluster(0);
        cluster1.addVariant(var1);

        tester.addClusterAndSVs(cluster1);
        tester.preClusteringInit();
        FoldbackFinder.markFoldbacks(tester.Analyser.getState().getChrBreakendMap());

        assertEquals(var1.getFoldbackId(true), var1.id());
        assertEquals(var1.getFoldbackId(false), var1.id());
        assertEquals(var1.getFoldbackLength(true), 100);

        // now test with a deletion bridge obscuring the backmost breakend
        tester.clearClustersAndSVs();

        SvCluster cluster2 = new SvCluster(0);
        cluster2.addVariant(var1);

        SvVarData dup = createDup(tester.nextVarId(), "1", 190, 400);
        cluster2.addVariant(dup);

        tester.addClusterAndSVs(cluster2);
        tester.preClusteringInit();
        FoldbackFinder.markFoldbacks(tester.Analyser.getState().getChrBreakendMap());

        assertEquals(var1.getFoldbackId(true), var1.id());
        assertEquals(var1.getFoldbackId(false), var1.id());
        assertEquals(var1.getFoldbackLength(true), 100);

        // now invalidate the foldback by putting a deletion bridge at the start
        tester.clearClustersAndSVs();

        SvVarData var2 = createInv(tester.nextVarId(), "1", 100, 200, -1);

        SvCluster cluster3 = new SvCluster(1);
        cluster3.addVariant(var2);

        SvVarData del = createDel(tester.nextVarId(), "1", 190, 400);
        cluster3.addVariant(del);

        tester.addClusterAndSVs(cluster3);
        tester.preClusteringInit();

        SvLinkedPair dbPair = var2.getDBLink(false);
        assertTrue(dbPair != null);
        assertTrue(dbPair.hasBreakend(var2, false));
        assertTrue(dbPair.hasBreakend(del, true));

        FoldbackFinder.markFoldbacks(tester.Analyser.getState().getChrBreakendMap());

        assertEquals(var2.getFoldbackBreakend(true), null);
        assertEquals(var2.getFoldbackBreakend(false), null);

        // finally an INV with DBs at both ends, making it invalid
        tester.clearClustersAndSVs();

        SvVarData var3 = createInv(tester.nextVarId(), "1", 100, 200, 1);

        SvCluster cluster4 = new SvCluster(2);
        cluster4.addVariant(var3);

        SvVarData sgl = createSgl(tester.nextVarId(), "1", 90, -1);
        cluster4.addVariant(sgl);

        SvVarData sgl2 = createSgl(tester.nextVarId(), "1", 190, -1);
        cluster4.addVariant(sgl2);

        tester.addClusterAndSVs(cluster4);
        tester.preClusteringInit();

        dbPair = var3.getDBLink(true);
        assertTrue(dbPair != null);
        assertTrue(dbPair.hasBreakend(var3, true));
        assertTrue(dbPair.hasBreakend(sgl, true));

        dbPair = var3.getDBLink(false);
        assertTrue(dbPair != null);
        assertTrue(dbPair.hasBreakend(var3, false));
        assertTrue(dbPair.hasBreakend(sgl2, true));

        FoldbackFinder.markFoldbacks(tester.Analyser.getState().getChrBreakendMap());

        assertEquals(var3.getFoldbackBreakend(true), null);
        assertEquals(var3.getFoldbackBreakend(false), null);
    }

}