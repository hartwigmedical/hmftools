package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_FOLDBACK;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_FOLDBACK_DSB;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_ISOLATED_BE;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_TI_ONLY;
import static com.hartwig.hmftools.linx.types.ArmCluster.getArmClusterData;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.analysis.ClusterMetrics;
import com.hartwig.hmftools.linx.types.ArmCluster;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import com.hartwig.hmftools.linx.utils.LinxTester;

public class AnnotationTest
{
    @Test
    public void testLocalTopology()
    {
        LinxTester tester = new LinxTester();

        // create a set of SVs which exhibit the various types of local topology
        final SvVarData var1 = createInv(1, "1", 101000, 107000, 1);

        // remote TI from BND, with other ends forming a DB and an isolated BE
        final SvVarData var2 = createBnd(2, "1", 149000, 1, "2", 1000, -1);
        final SvVarData var3 = createBnd(3, "1", 250000, 1, "2", 1100, 1);

        // DUP forms a DB with the next foldback and a DB with the previous BND
        final SvVarData var4 = createDup(4, "1", 150000, 200500);
        final SvVarData var5 = createInv(5, "1", 201000, 204000, -1);

        // complex foldback
        final SvVarData var6 = createInv(6, "1", 301000, 304000, 1);
        final SvVarData var7 = createDel(7, "1", 305000, 306000);
        final SvVarData var8 = createInv(8, "1", 310000, 315000, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);
        tester.AllVariants.add(var8);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.isFoldback());
        assertTrue(var5.isFoldback());
        assertTrue(var6.isFoldback());
        assertTrue(var8.isFoldback());

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        final int[] armClusterData = getArmClusterData(cluster);

        assertEquals(6, cluster.getArmClusters().size());
        assertEquals(1, armClusterData[ARM_CL_ISOLATED_BE]);
        assertEquals(1, armClusterData[ARM_CL_TI_ONLY]);
        assertEquals(1, armClusterData[ARM_CL_DSB]);
        assertEquals(1, armClusterData[ARM_CL_FOLDBACK]);
        assertEquals(1, armClusterData[ARM_CL_FOLDBACK_DSB]);
        assertEquals(1, armClusterData[ARM_CL_COMPLEX_FOLDBACK]);
    }

    @Test
    public void testMultipleTIs()
    {
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 1000, -1);
        final SvVarData var2 = createBnd(tester.nextVarId(), "1", 1110, -1, "2", 1200, 1);
        final SvVarData var3 = createBnd(tester.nextVarId(), "1", 1310, 1, "2", 1210, -1);
        final SvVarData var4 = createBnd(tester.nextVarId(), "1", 1320, -1, "2", 1410, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        final int[] armClusterData = getArmClusterData(cluster);

        assertEquals(2, cluster.getArmClusters().size());
        assertEquals(1, armClusterData[ARM_CL_TI_ONLY]);
        assertEquals(1, armClusterData[ARM_CL_DSB]);

        long armClusterTIs = cluster.getArmClusters().stream().mapToInt(x -> x.getTICount()).sum();
        assertEquals(3, armClusterTIs);

        boolean hasDoubleTI = false;

        for(ArmCluster armCluster : cluster.getArmClusters())
        {
            if(armCluster.getTICount() == 2 && armCluster.getType() == ARM_CL_TI_ONLY)
                hasDoubleTI = true;
        }

        assertTrue(hasDoubleTI);
    }

    @Test
    public void testClusterMetrics()
    {
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createInv(tester.nextVarId(), "1", 1000, 49000,1);
        final SvVarData var2 = createInv(tester.nextVarId(), "1", 10000, 50000,-1);

        // DELs unclustered over TI
        final SvVarData del1 = createDel(tester.nextVarId(), "1", 20000, 20200);
        final SvVarData del2 = createDel(tester.nextVarId(), "1", 30000, 30200);

        final SvVarData bnd1 = createBnd(tester.nextVarId(), "1", 52000, 1, "2", 100, -1);
        final SvVarData bnd2 = createBnd(tester.nextVarId(), "1", 52050, -1, "2", 200, 1); // short DB

        final SvVarData del3 = createDel(tester.nextVarId(), "1", 60000, 60100); // unclustered over TI

        final SvVarData bnd3 = createBnd(tester.nextVarId(), "1", 74000, 1, "2", 10000, -1);
        final SvVarData bnd4 = createBnd(tester.nextVarId(), "1", 75000, -1, "2", 11000, 1);

        final SvVarData var3 = createDel(tester.nextVarId(), "1", 76000, 90000);
        final SvVarData var4 = createDup(tester.nextVarId(), "1", 77000, 90100);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(del1);
        tester.AllVariants.add(del2);
        tester.AllVariants.add(del3);
        tester.AllVariants.add(bnd1);
        tester.AllVariants.add(bnd2);
        tester.AllVariants.add(bnd3);
        tester.AllVariants.add(bnd4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(4, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        final ClusterMetrics metrics = cluster.getMetrics();
        metrics.populate(cluster, tester.Analyser.getState().getChrBreakendMap());

        assertEquals(5, metrics.DBCount);
        assertEquals(1, metrics.ShortDBCount);
        assertEquals(76000, metrics.TotalRange);
        assertEquals(3, metrics.TraversedDelCount);
        assertEquals(1, metrics.TraversedShortDelCount);
        assertEquals(5, metrics.ImpliedTICount);
        assertEquals(500, metrics.TraversedDelLength);
    }

}
