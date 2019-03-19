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

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLOH;
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

        SvVarData var6 = createBnd(tester.nextVarId(), "1", 6000, -1, "5", 1000, 1);
        tester.AllVariants.add(var6);

        // equivalent breakends are kept separate
        SvVarData var7 = createSgl(tester.nextVarId(), "1", 7000, -1, false);
        tester.AllVariants.add(var7);

        SvVarData var8 = createSgl(tester.nextVarId(), "1", 7000, -1, false);
        tester.AllVariants.add(var8);

        SvVarData var9 = createSgl(tester.nextVarId(), "1", 7100, -1, false);
        var9.setAssemblyData(true, ASSEMBLY_TYPE_EQV);
        tester.AllVariants.add(var9);

        // and some variants on the BND's other chromosome, and other linking BNDs
        SvVarData var10 = createIns(tester.nextVarId(), "5", 2000, 20000);
        tester.AllVariants.add(var10);

        SvVarData var11 = createBnd(tester.nextVarId(), "2", 6000, -1, "5", 21000, 1);
        tester.AllVariants.add(var11);

        SvVarData var12 = createIns(tester.nextVarId(), "2", 7000, 8000);
        tester.AllVariants.add(var12);

        SvVarData var13 = createBnd(tester.nextVarId(), "3", 6000, -1, "5", 22000, 1);
        tester.AllVariants.add(var13);

        tester.preClusteringInit();

        tester.mergeOnProximity();

        assertEquals(tester.getClusters().size(), 4);
        assertEquals(tester.getClusters().get(0).getSvCount(), 10);
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

        tester.preClusteringInit();

        tester.mergeOnProximity();

        assertEquals(tester.getClusters().size(), 2);
        assertTrue(tester.getClusters().get(0).getSVs().contains(del1));
        assertTrue(tester.getClusters().get(0).getSVs().contains(del3));
        assertTrue(tester.getClusters().get(1).getSVs().contains(del2));
        assertTrue(tester.getClusters().get(1).getSVs().contains(del4));
    }

    @Test
    public void testFoldbackMerge()
    {
        SvTestHelper tester = new SvTestHelper();

        // 2 clusters with foldbacks on the same arm are merged
        SvVarData inv1 = createInv(tester.nextVarId(), "1", 100, 200, -1);
        tester.AllVariants.add(inv1);

        SvVarData inv2 = createInv(tester.nextVarId(), "1", 20000, 20100, 1);
        tester.AllVariants.add(inv2);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.getClusters().size());
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

        tester.Analyser.clusterAndAnalyse();

        assertEquals(tester.getClusters().size(), 2);
        assertTrue(tester.getClusters().get(0).getSVs().contains(inv2));
        assertTrue(tester.getClusters().get(0).getSVs().contains(sgl));
        assertTrue(tester.getClusters().get(1).getSVs().contains(del));

        assertTrue(inv1.getClusterReason().contains(CLUSTER_REASON_FOLDBACKS));
        assertTrue(inv2.getClusterReason().contains(CLUSTER_REASON_FOLDBACKS));

        tester.clearClustersAndSVs();

        // test with a single foldback facing the centromere
        SvVarData inv3 = createInv(tester.nextVarId(), "2", 10000, 10100, -1);
        tester.AllVariants.add(inv3);

        SvVarData sgl2 = createSgl(tester.nextVarId(), "2", 20000, 1, false);
        tester.AllVariants.add(sgl2);

        // next is too far away
        SvVarData sgl3 = createSgl(tester.nextVarId(), "2", 10000000, 1, false);
        tester.AllVariants.add(sgl3);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(tester.getClusters().size(), 2);
        assertTrue(tester.getClusters().get(0).getSVs().contains(inv3));
        assertTrue(tester.getClusters().get(0).getSVs().contains(sgl2));
        assertTrue(tester.getClusters().get(1).getSVs().contains(sgl3));

        assertTrue(inv3.getClusterReason().contains(CLUSTER_REASON_FOLDBACKS));
        assertTrue(sgl2.getClusterReason().contains(CLUSTER_REASON_FOLDBACKS));

    }

    @Test
    public void testLohResolvingClusterMerge()
    {
        SvTestHelper tester = new SvTestHelper();

        SvVarData var1 = createDup("1", "1", 10000, 50000);
        SvVarData var2 = createBnd("2", "1", 20000, 1, "3", 1000, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        Map<String, List<SvLOH>> lohDataMap = new HashMap();
        List<SvLOH> lohData = Lists.newArrayList();

        lohData.add(new SvLOH(tester.SampleId, "1", 1, 2, 50000, 1000000,
                "DUP", "TELOMERE", 1, 1, 1, 0, 1, 1,
                var1.id(), "", false, true));

        lohDataMap.put(tester.SampleId, lohData);

        tester.ClusteringMethods.setSampleLohData(lohDataMap);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // now check final chain-finding across all sub-clusters
        assertEquals(tester.Analyser.getClusters().size(), 1);
    }

    @Test
    public void testConsistentBreakendOverlapMerge()
    {
        SvTestHelper tester = new SvTestHelper();

        List<SvVarData> allVariants = Lists.newArrayList();

        // a cluster has 3 consecutive breakends which span other unresolved SV breakends, which are then merged in
        SvVarData consec1 = createBnd(tester.nextVarId(), "1", 100000, 1, "2", 100, -1);
        allVariants.add(consec1);

        SvVarData consec2 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 200, 1);
        allVariants.add(consec2);

        SvVarData consec3 = createInv(tester.nextVarId(), "1", 30000, 101000, 1);
        allVariants.add(consec3);

        // create some SV in resolved clusters which will be ignored - a DEL with external TI , a simple DEL and a low-qual
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 10000, 1, "3", 200, 1);
        allVariants.add(var1);

        SvVarData var2 = createBnd(tester.nextVarId(), "1", 10100, -1, "3", 100, -1);
        allVariants.add(var2);

        SvVarData var3 = createDel(tester.nextVarId(), "1", 60000, 60100);
        allVariants.add(var3);

        SvVarData var4 = createSgl(tester.nextVarId(), "1", 80000, -1, false);
        var4.setAssemblyData(true, ASSEMBLY_TYPE_EQV);
        allVariants.add(var4);

        // now some SV which will be overlapped by the consecutive breakends
        SvVarData overlap1 = createBnd(tester.nextVarId(), "1", 20000, -1, "4", 200, 1);
        allVariants.add(overlap1);

        SvVarData overlap2 = createSgl(tester.nextVarId(), "1", 40000, -1, false);
        allVariants.add(overlap2);

        allVariants.addAll(allVariants);

        tester.AllVariants.addAll(allVariants);
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(4, tester.getClusters().size());

        SvCluster mainCluster = null;
        for(final SvCluster cluster : tester.getClusters())
        {
            if(cluster.getSvCount() == 5)
            {
                mainCluster = cluster;
                break;
            }
        }

        if(mainCluster == null)
            assertTrue(false);

        assertTrue(consec1.getClusterReason().contains(CLUSTER_REASON_LOOSE_OVERLAP));
        assertTrue(overlap1.getClusterReason().contains(CLUSTER_REASON_LOOSE_OVERLAP));
        assertTrue(overlap2.getClusterReason().contains(CLUSTER_REASON_LOOSE_OVERLAP));
    }

}
