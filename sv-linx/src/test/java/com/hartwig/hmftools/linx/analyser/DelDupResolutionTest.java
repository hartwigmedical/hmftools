package com.hartwig.hmftools.linx.analyser;


import static com.hartwig.hmftools.linx.types.ResolvedType.DUP_TI;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import static com.hartwig.hmftools.linx.utils.SvTestRoutines.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestRoutines.createDup;

import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import com.hartwig.hmftools.linx.utils.LinxTester;

public class DelDupResolutionTest
{

    @Test
    public void testEnclosedDelInDup()
    {
        LinxTester tester = new LinxTester();

        SvVarData var1 = createDup(tester.nextVarId(), "1", 1000, 50000);
        SvVarData var2 = createDel(tester.nextVarId(), "1", 5000, 40000);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(cluster.getResolvedType() == DUP_TI);

        tester.clearClustersAndSVs();

        /*

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
        assertEquals(cluster.getSyntheticLength(), var1.position(true) - var4.position(true));
        assertEquals(600, cluster.getSyntheticTILength()); // takes average

        tester.clearClustersAndSVs();

        // special case of the reciprocal DUP which is more likely to not be chained through a long TI, and so is split
        var1 = createBnd(tester.nextVarId(), "1", 20000, 1, "3", 1000, -1);
        var2 = createBnd(tester.nextVarId(), "2", 1000, -1, "3", 1500, 1);

        var3 = createBnd(tester.nextVarId(), "1", 1000, -1, "3", 2000, -1);
        var4 = createBnd(tester.nextVarId(), "2", 20000, 1, "3", 2500, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);
        assertEquals(2, cluster.getChains().size());

        assertTrue(!cluster.isResolved());
        assertEquals(RECIP_DUPS, cluster.getResolvedType());
        assertTrue(cluster.isSyntheticType());
        assertEquals(0, cluster.getSyntheticLength());
        assertEquals(500, cluster.getSyntheticTILength());

        // test the 2 INV (typically foldbacks) case

        tester.clearClustersAndSVs();

        var1 = createInv(tester.nextVarId(), "1", 1000, 3000, -1);
        var2 = createInv(tester.nextVarId(), "1", 51000, 53000, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);
        assertEquals(1, cluster.getChains().size());

        assertTrue(!cluster.isResolved());
        assertEquals(FB_INV_PAIR, cluster.getResolvedType());
        assertTrue(!cluster.isSyntheticType());
        assertEquals(var2.position(false) - var1.position(true), cluster.getSyntheticLength());
        assertEquals(var2.position(true) - var1.position(false), cluster.getSyntheticTILength());

        // and finally a local version of the overlapping DUPs scenarios with a reconfigured chain
        tester.clearClustersAndSVs();

        var1 = createInv(tester.nextVarId(), "1", 1000, 51000, -1);
        var2 = createInv(tester.nextVarId(), "1", 5000, 53000, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);
        assertEquals(1, cluster.getChains().size());

        assertTrue(!cluster.isResolved());
        assertEquals(RECIP_DUPS, cluster.getResolvedType());
        assertTrue(!cluster.isSyntheticType());
        assertEquals(var1.position(false) - var2.position(true), cluster.getSyntheticLength());
        assertEquals(var2.position(false) - var1.position(true), cluster.getSyntheticTILength());

        // and test again with additional short TIs in the mix
        tester.clearClustersAndSVs();

        var1 = createBnd(tester.nextVarId(), "1", 1000, -1, "2", 2000, -1);
        var2 = createBnd(tester.nextVarId(), "1", 51000, -1, "2", 2500, 1);

        var3 = createBnd(tester.nextVarId(), "1", 5000, 1, "2", 4000, -1);
        var4 = createBnd(tester.nextVarId(), "1", 53000, 1, "2", 4500, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);
        assertEquals(1, cluster.getChains().size());

        assertTrue(!cluster.isResolved());
        assertEquals(RECIP_DUPS, cluster.getResolvedType());
        assertTrue(cluster.isSyntheticType());
        assertEquals(var2.position(true) - var3.position(true), cluster.getSyntheticLength());
        assertEquals(var4.position(true) - var1.position(true), cluster.getSyntheticTILength());

        */
    }

}
