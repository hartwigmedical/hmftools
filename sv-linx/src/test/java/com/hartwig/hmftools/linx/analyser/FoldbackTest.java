package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSgl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.analysis.FoldbackFinder;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class FoldbackTest
{

    @Test
    public void testInvFoldbacks()
    {
        LinxTester tester = new LinxTester();

        // create foldbacks from a simple inversion

        SvVarData var1 = createInv(tester.nextVarId(), "1", 100, 200, 1);

        SvCluster cluster = new SvCluster(0);
        cluster.addVariant(var1);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();
        FoldbackFinder.markFoldbacks(tester.Analyser.getState().getChrBreakendMap());

        assertEquals(var1.getFoldbackId(true), var1.id());
        assertEquals(var1.getFoldbackId(false), var1.id());
        assertEquals(var1.getFoldbackLength(true), 100);

        // now test with a deletion bridge obscuring the backmost breakend
        tester.clearClustersAndSVs();

        cluster = new SvCluster(0);
        cluster.addVariant(var1);

        SvVarData dup = createDup(tester.nextVarId(), "1", 190, 400);
        cluster.addVariant(dup);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();
        FoldbackFinder.markFoldbacks(tester.Analyser.getState().getChrBreakendMap());

        assertEquals(var1.getFoldbackId(true), var1.id());
        assertEquals(var1.getFoldbackId(false), var1.id());
        assertEquals(var1.getFoldbackLength(true), 100);

        // and again with the foldback facing the other way
        tester.clearClustersAndSVs();

        cluster = new SvCluster(0);

        SvVarData var2 = createInv(tester.nextVarId(), "1", 100, 200, -1);
        cluster.addVariant(var2);

        SvVarData sgl = createSgl(tester.nextVarId(), "1", 110, 1);
        cluster.addVariant(sgl);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();
        FoldbackFinder.markFoldbacks(tester.Analyser.getState().getChrBreakendMap());

        assertEquals(var2.getFoldbackId(true), var2.id());
        assertEquals(var2.getFoldbackId(false), var2.id());
        assertEquals(var2.getFoldbackLength(true), 100);

        // now invalidate the foldback by putting a deletion bridge at the start
        tester.clearClustersAndSVs();

        var2 = createInv(tester.nextVarId(), "1", 100, 200, -1);

        cluster = new SvCluster(1);
        cluster.addVariant(var2);

        SvVarData del = createDel(tester.nextVarId(), "1", 190, 400);
        cluster.addVariant(del);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();

        DbPair dbPair = var2.getDBLink(false);
        assertTrue(dbPair != null);
        assertTrue(dbPair.hasBreakend(var2, false));
        assertTrue(dbPair.hasBreakend(del, true));

        FoldbackFinder.markFoldbacks(tester.Analyser.getState().getChrBreakendMap());

        assertEquals(var2.getFoldbackBreakend(true), null);
        assertEquals(var2.getFoldbackBreakend(false), null);

        // finally an INV with DBs at both ends, making it invalid
        tester.clearClustersAndSVs();

        SvVarData var3 = createInv(tester.nextVarId(), "1", 100, 200, 1);

        cluster = new SvCluster(2);
        cluster.addVariant(var3);

        sgl = createSgl(tester.nextVarId(), "1", 90, -1);
        cluster.addVariant(sgl);

        SvVarData sgl2 = createSgl(tester.nextVarId(), "1", 190, -1);
        cluster.addVariant(sgl2);

        tester.addClusterAndSVs(cluster);
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

        assertEquals(null, var3.getFoldbackBreakend(true));
        assertEquals(null, var3.getFoldbackBreakend(false));
    }

}