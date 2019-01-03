package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createBnd;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createIns;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createSgl;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_FOLDBACKS;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_LOOSE_OVERLAP;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.ASSEMBLY_TYPE_EQV;

import static org.junit.Assert.assertEquals;

import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Ignore;
import org.junit.Test;

public class MergeRuleTests
{

    @Test
    public void testProximityMerge()
    {
        SvTestHelper tester = new SvTestHelper();

        // basic proximity clustering
        SvVarData var1 = createDel(tester.nextVarId(), "1", 1000, 1100);
        tester.AllVariants.add(var1);

        SvVarData var2 = createDup(tester.nextVarId(), "1", 2000, 2100);
        tester.AllVariants.add(var2);

        SvVarData var3 = createIns(tester.nextVarId(), "1", 3000, 3100);
        tester.AllVariants.add(var3);

        SvVarData var4 = createInv(tester.nextVarId(), "1", 4000, 4100, 1);
        tester.AllVariants.add(var4);

        SvVarData var5 = createSgl(tester.nextVarId(), "1", 5000, -1, false);
        tester.AllVariants.add(var5);

        SvVarData var6 = createBnd(tester.nextVarId(), "1", 6000, -1, "2", 1000, 1);
        tester.AllVariants.add(var6);

        // equivalent breakends are kept separate
        SvVarData var7 = createSgl(tester.nextVarId(), "1", 7000, -1, false);
        tester.AllVariants.add(var7);

        SvVarData var8 = createSgl(tester.nextVarId(), "1", 7000, -1, false);
        tester.AllVariants.add(var8);

        SvVarData var9 = createSgl(tester.nextVarId(), "1", 7100, -1, false);
        var9.setAssemblyData(true, ASSEMBLY_TYPE_EQV);
        tester.AllVariants.add(var9);

        tester.preClusteringInit();

        tester.mergeOnProximity();

        assertEquals(tester.getClusters().size(), 4);
        assertEquals(tester.getClusters().get(0).getCount(), 6);
        assertTrue(tester.getClusters().get(1).getSVs().contains(var7));
        assertTrue(tester.getClusters().get(2).getSVs().contains(var8));
        assertTrue(tester.getClusters().get(3).getSVs().contains(var9));

        // non-overlapping DELs can merge, whereas overlapping DELs are split out
        tester.clearClustersAndSVs();

        SvVarData del1 = createDel(tester.nextVarId(), "2", 1000, 3000);
        tester.AllVariants.add(del1);

        SvVarData del2 = createDel(tester.nextVarId(), "2", 2000, 4000);
        tester.AllVariants.add(del2);

        SvVarData del3 = createDel(tester.nextVarId(), "2", 3500, 8000);
        tester.AllVariants.add(del3);

        SvVarData del4 = createDel(tester.nextVarId(), "2", 6000, 10000);
        tester.AllVariants.add(del4);

        tester.mergeOnProximity();

        assertEquals(tester.getClusters().size(), 2);
        assertTrue(tester.getClusters().get(0).getSVs().contains(del1));
        assertTrue(tester.getClusters().get(0).getSVs().contains(del3));
        assertTrue(tester.getClusters().get(0).getSVs().contains(del4));
        assertTrue(tester.getClusters().get(1).getSVs().contains(del2));

        // test exclusion of low-quality SVs
        tester.clearClustersAndSVs();
    }

    @Test
    @Ignore
    public void testFoldbackMerge()
    {
        SvTestHelper tester = new SvTestHelper();

        // 2 clusters with foldbacks on the same arm are merged
        SvVarData inv1 = createInv(tester.nextVarId(), "1", 100, 200, -1);
        tester.AllVariants.add(inv1);

        SvVarData inv2 = createInv(tester.nextVarId(), "1", 20000, 20100, 1);
        tester.AllVariants.add(inv2);

        tester.preClusteringInit();

        tester.mergeOnProximity();
        assertEquals(tester.getClusters().size(), 2);

        tester.Analyser.clusterAndAnalyse();

        assertEquals(tester.getClusters().size(), 1);
        assertTrue(inv1.getClusterReason().contains(CLUSTER_REASON_FOLDBACKS));
        assertTrue(inv2.getClusterReason().contains(CLUSTER_REASON_FOLDBACKS));

        // non-overlapping DELs can merge, whereas overlapping DELs are split out
        tester.clearClustersAndSVs();

        // test again with a foldback and an opposing SV as the next breakend
        tester.AllVariants.add(inv2);

        // the DEL is resolved and will be ignored
        SvVarData del = createDel(tester.nextVarId(), "1", 10000, 10100);
        tester.AllVariants.add(del);

        SvVarData sgl = createSgl(tester.nextVarId(), "1", 1000, -1, false);
        tester.AllVariants.add(sgl);

        tester.preClusteringInit();

        tester.mergeOnProximity();
        assertEquals(tester.getClusters().size(), 3);

        tester.Analyser.clusterAndAnalyse();

        assertEquals(tester.getClusters().size(), 2);
        assertTrue(tester.getClusters().get(0).getSVs().contains(inv2));
        assertTrue(tester.getClusters().get(0).getSVs().contains(sgl));
        assertTrue(tester.getClusters().get(1).getSVs().contains(del));

        assertTrue(inv1.getClusterReason().contains(CLUSTER_REASON_FOLDBACKS));
        assertTrue(inv2.getClusterReason().contains(CLUSTER_REASON_FOLDBACKS));

    }

    @Test
    public void testConsistentBreakendOverlapMerge()
    {
        SvTestHelper tester = new SvTestHelper();

        // a cluster has 3 consecutive breakends which span other unresolved SV breakends, which are then merged in
        SvVarData consec1 = createBnd(tester.nextVarId(), "1", 100000, 1, "2", 100, -1);
        tester.AllVariants.add(consec1);

        SvVarData consec2 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 200, 1);
        tester.AllVariants.add(consec2);

        SvVarData consec3 = createInv(tester.nextVarId(), "1", 30000, 101000, 1);
        tester.AllVariants.add(consec3);

        // create some SV in resolved clusters which will be ignored - a DEL with external TI , a simple DEL and a low-qual
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 10000, 1, "3", 200, 1);
        tester.AllVariants.add(var1);

        SvVarData var2 = createBnd(tester.nextVarId(), "1", 10100, -1, "3", 100, -1);
        tester.AllVariants.add(var2);

        SvVarData var3 = createDel(tester.nextVarId(), "1", 60000, 60100);
        tester.AllVariants.add(var3);

        SvVarData var4 = createSgl(tester.nextVarId(), "1", 80000, -1, false);
        var4.setAssemblyData(true, ASSEMBLY_TYPE_EQV);
        tester.AllVariants.add(var4);

        // now some SV which will be overlapped by the consecutive breakends
        SvVarData overlap1 = createBnd(tester.nextVarId(), "1", 20000, -1, "4", 200, 1);
        tester.AllVariants.add(overlap1);

        SvVarData overlap2 = createSgl(tester.nextVarId(), "1", 40000, -1, false);
        tester.AllVariants.add(overlap2);

        tester.preClusteringInit();

        tester.mergeOnProximity();
        assertEquals(tester.getClusters().size(), 6);

        tester.Analyser.clusterAndAnalyse();

        assertEquals(tester.getClusters().size(), 4);
        final SvCluster mainCluster = tester.getClusters().get(0);
        assertEquals(mainCluster.getCount(), 5);
        assertTrue(consec1.getClusterReason().contains(CLUSTER_REASON_LOOSE_OVERLAP));
        assertTrue(overlap1.getClusterReason().contains(CLUSTER_REASON_LOOSE_OVERLAP));
        assertTrue(overlap2.getClusterReason().contains(CLUSTER_REASON_LOOSE_OVERLAP));
    }

}
