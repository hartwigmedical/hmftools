package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.linx.analysis.ClusterClassification.getSyntheticLength;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.getSyntheticTiLength;
import static com.hartwig.hmftools.linx.types.ArmClusterType.COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.linx.types.ArmClusterType.COMPLEX_OTHER;
import static com.hartwig.hmftools.linx.types.ArmClusterType.DSB;
import static com.hartwig.hmftools.linx.types.ArmClusterType.TI_ONLY;
import static com.hartwig.hmftools.linx.types.ResolvedType.DEL;
import static com.hartwig.hmftools.linx.types.ResolvedType.DUP;
import static com.hartwig.hmftools.linx.types.ResolvedType.FB_INV_PAIR;
import static com.hartwig.hmftools.linx.types.ResolvedType.INV;
import static com.hartwig.hmftools.linx.types.ResolvedType.PAIR_OTHER;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV;
import static com.hartwig.hmftools.linx.types.ResolvedType.SGL;
import static com.hartwig.hmftools.linx.types.ResolvedType.UNBAL_TRANS;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSgl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class SyntheticTest
{
    @Test
    public void testSyntheticDelDupFromInvPairs()
    {
        LinxTester tester = new LinxTester();

        // create 2 INVs with varying positions to check what synthetic DEL or DUP they create

        // no overlap but 2 deletion bridges - doesn't make anything
        SvVarData var1 = createInv(tester.nextVarId(), "1", 100, 200, 1);
        SvVarData var2 = createInv(tester.nextVarId(), "1", 300, 400, -1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(!cluster.isResolved());
        assertTrue(cluster.getResolvedType() == PAIR_OTHER);

        // 1 pair of overlapping breakends
        var1 = createInv(tester.nextVarId(), "1", 100, 300, 1);
        var2 = createInv(tester.nextVarId(), "1", 200, 400, -1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DEL); // since TI is short

        assertEquals(1, cluster.getArmClusters().size());
        assertEquals(DSB, cluster.getArmClusters().get(0).getType());
        assertEquals(1, cluster.getArmClusters().get(0).getTICount());

        // test 2 DSBs but with an overlapping end less than permitted TI length
        var1 = createInv(tester.nextVarId(), "1", 100, 2400, 1);
        var2 = createInv(tester.nextVarId(), "1", 90, 3000, -1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RECIP_INV);

        // test 2 overlapping breakends
        var1 = createInv(tester.nextVarId(), "1", 100, 10400, 1);
        var2 = createInv(tester.nextVarId(), "1", 250, 10350, -1);
        var1.setAssemblyData(false, "asmb12");
        var2.setAssemblyData(false, "asmb12");

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DEL);
        assertEquals(getSyntheticLength(cluster), var2.position(true) - var1.position(true));
        assertEquals(getSyntheticTiLength(cluster), var1.position(false) - var2.position(false) + 1);

        assertEquals(2, cluster.getArmClusters().size());
        assertEquals(DSB, cluster.getArmClusters().get(0).getType());
        assertEquals(TI_ONLY, cluster.getArmClusters().get(1).getType());

        // test 2 overlapping breakends but where a pair of breakend form an overlapping DB
        var1 = createInv(tester.nextVarId(), "1", 500, 4000, 1);
        var2 = createInv(tester.nextVarId(), "1", 480, 5000, -1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RECIP_INV);

        // 3 overlapping breakends
        var1 = createInv(tester.nextVarId(), "1", 200, 10400, 1);
        var2 = createInv(tester.nextVarId(), "1", 100, 10350, -1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DUP);
        assertEquals(getSyntheticLength(cluster), var1.position(true) - var2.position(true));
        assertEquals(getSyntheticTiLength(cluster), var1.position(false) - var2.position(false) + 1);

        assertEquals(2, cluster.getArmClusters().size());
        assertEquals(COMPLEX_OTHER, cluster.getArmClusters().get(0).getType());
        assertEquals(TI_ONLY, cluster.getArmClusters().get(1).getType());

        // 4 overlapping breakends - also forming 2 facing foldbacks
        var1 = createInv(tester.nextVarId(), "1", 2300, 2400, 1);
        var2 = createInv(tester.nextVarId(), "1", 100, 250, -1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(!cluster.isResolved());
        assertTrue(cluster.getResolvedType() == FB_INV_PAIR);

        assertEquals(1, cluster.getArmClusters().size());
        assertEquals(COMPLEX_FOLDBACK, cluster.getArmClusters().get(0).getType());
    }

    @Test
    public void testSyntheticDelDupFromBndPairs()
    {
        LinxTester tester = new LinxTester();

        // create 2 BNDs with varying positions to check what synthetic DEL or DUP they create

        // single TI with a DEL
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 100, 1, "2", 200, 1);
        SvVarData var2 = createBnd(tester.nextVarId(), "1", 200, -1, "2", 150, -1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DEL);
        assertEquals(getSyntheticLength(cluster), var2.position(true) - var1.position(true));
        assertEquals(getSyntheticTiLength(cluster), var1.position(false) - var2.position(false) + 1);

        // 2 TIs with a DUP
        var1 = createBnd(tester.nextVarId(), "1", 200, 1, "2", 200, 1);
        var2 = createBnd(tester.nextVarId(), "1", 100, -1, "2", 150, -1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DUP);
        assertEquals(getSyntheticLength(cluster), var1.position(true) - var2.position(true));
        assertEquals(getSyntheticTiLength(cluster), var1.position(false) - var2.position(false) + 1);
    }

    @Test
    public void testSyntheticDelDupFromDelsAndDups()
    {
        LinxTester tester = new LinxTester();

        // DEL and DUP - DEL enclosed
        SvVarData var1 = createDel(tester.nextVarId(), "1", 300, 350);
        SvVarData var2 = createDup(tester.nextVarId(), "1", 100, 400);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DUP);
        assertEquals(getSyntheticLength(cluster), var1.position(true) - var2.position(true));
        assertEquals(getSyntheticTiLength(cluster), var2.position(false) - var1.position(false) + 1);

        // DEL and DUP - DEL overlapping DUP
        var1 = createDel(tester.nextVarId(), "1", 100, 350);
        var2 = createDup(tester.nextVarId(), "1", 200, 400);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DEL);
        assertEquals(getSyntheticLength(cluster), var2.position(true) - var1.position(true));
        assertEquals(getSyntheticTiLength(cluster), var2.position(false) - var1.position(false) + 1);

        // DEL and DUP - overlapping the other way
        var1 = createDup(tester.nextVarId(), "1", 100, 350);
        var2 = createDel(tester.nextVarId(), "1", 200, 400);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DEL);
        assertEquals(getSyntheticLength(cluster), var2.position(false) - var1.position(false));
        assertEquals(getSyntheticTiLength(cluster), var2.position(true) - var1.position(true) + 1);

        // 2 DUPs next to each other - since the short TI overlaps the chain ends, it's not considered a synthetic
        var1 = createDup(tester.nextVarId(), "1", 100, 2200);
        var2 = createDup(tester.nextVarId(), "1", 300, 2400);

        var1.getTIAssemblies(true).add("asmb1");
        var2.getTIAssemblies(false).add("asmb1");

        tester.addAndCluster(var1, var2);

        assertEquals(2, tester.Analyser.getClusters().size());

        // 2 DUPs with an overlap
        var1 = createDup(tester.nextVarId(), "1", 100, 300);
        var2 = createDup(tester.nextVarId(), "1", 250, 400);

        var1.getTIAssemblies(false).add("asmb1");
        var2.getTIAssemblies(true).add("asmb1");

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DUP);
    }

    @Test
    public void testSyntheticDelsAndDupsWithLongerChains()
    {
        // as above but with more than 1 short TI
        LinxTester tester = new LinxTester();

        SvVarData var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 1000, -1);
        SvVarData var2 = createBnd(tester.nextVarId(), "2", 1200, 1, "3", 2000, 1);
        SvVarData var3 = createBnd(tester.nextVarId(), "1", 5000, -1, "3", 1600, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DEL);
        assertEquals(getSyntheticLength(cluster), var3.position(true) - var1.position(true));
        assertEquals(401, getSyntheticTiLength(cluster)); // takes the longest

        tester.clearClustersAndSVs();

        // add one more SV to extend the chain into a DUP - the INV comes in too close for a TI with the first SV
        SvVarData var4 = createDup(tester.nextVarId(), "1", 980, 5600);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DUP);
        assertTrue(cluster.isSyntheticType());
        assertEquals(getSyntheticLength(cluster), var1.position(true) - var4.position(true));
        assertEquals(601, getSyntheticTiLength(cluster)); // takes average

        tester.clearClustersAndSVs();
    }

    @Test
    public void testSyntheticIncompletes()
    {
        LinxTester tester = new LinxTester();

        SvVarData var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 1000, -1);
        SvVarData var2 = createBnd(tester.nextVarId(), "2", 1200, 1, "3", 2000, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertFalse(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == UNBAL_TRANS);

        tester.clearClustersAndSVs();

        var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 1000, -1);
        var2 = createBnd(tester.nextVarId(), "1", 1200, 1, "2", 1500, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertFalse(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == INV);

        tester.clearClustersAndSVs();

        var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 1000, -1);
        var2 = createSgl(tester.nextVarId(), "2", 1500, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertFalse(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == SGL);
    }



}
