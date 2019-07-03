package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createSgl;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createTestSv;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.ALL_ANNOTATIONS;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNONTATION_DM;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

public class DoubleMinuteTest
{

    @Test
    public void testSimpleDupDM()
    {
        SvTestHelper tester = new SvTestHelper();

        // tester.logVerbose(true);
        tester.Config.RequiredAnnotations = ALL_ANNOTATIONS;

        // first a simple DUP
        final SvVarData dup1 = createTestSv("0","1","1",500,600,-1,1, DUP,10);
        dup1.setPloidyRecalcData(8, 12);

        tester.AllVariants.add(dup1);

        tester.addCopyNumberData();

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());

        SvCluster cluster = tester.Analyser.getClusters().get(0);
        assertEquals(1, cluster.getSvCount());
        assertTrue(cluster.getSVs().contains(dup1));
        assertTrue(cluster.getAnnotations().contains(CLUSTER_ANNONTATION_DM));
    }

    @Test
    public void testChainedDM()
    {
        // form a DM from 3 chained SVs, with some other SVs in the cluster having a lower ploidy
        SvTestHelper tester = new SvTestHelper();

        tester.logVerbose(true);
        tester.Config.RequiredAnnotations = ALL_ANNOTATIONS;

        // 1 s100 -> 6 e600-500s -> 4 s1500-100e -> 3 s200-2000e -> 5 s2100-2200e -> 2 e2500-1200s -> 1 e1000

        final SvVarData var1 = createTestSv("1","1","1",10100,11000,-1,-1, INV,10);
        final SvVarData var2 = createTestSv("2","1","1",11200,12500,1,1, INV,10);
        final SvVarData var3 = createTestSv("3","1","2",12000,10200,-1,1, BND,10);
        final SvVarData var4 = createTestSv("4","1","2",11500,10100,1,-1, BND,10);
        final SvVarData var5 = createTestSv("5","1","1",12100,12200,1,1, DEL,1);
        final SvVarData var6 = createTestSv("6","1","1",10500,10600,-1,1, DUP,2);

        // unrelated SVs at either end of the cluster
        final SvVarData other1 = createTestSv("7","1","1",100,200,1,-1, DEL,1);
        final SvVarData other2 = createTestSv("8","1","1",20000,20100,1,-1, DEL,1);
        final SvVarData other3 = createTestSv("9","2","2",20100,20100,1,-1, DEL,1);

        tester.AllVariants.add(other1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);

        tester.AllVariants.add(other2);
        tester.AllVariants.add(other3);
        tester.preClusteringInit();

        tester.addCopyNumberData();

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        boolean matched = false;

        for(final SvCluster cluster : tester.Analyser.getClusters())
        {
            if (cluster.getSvCount() == 6 && cluster.getAnnotations().contains(CLUSTER_ANNONTATION_DM))
            {
                matched = true;
                break;
            }
        }

        assertTrue(matched);
    }

    @Test
    public void testInvalidDM()
    {
        // first a cluster which grows in ploidy evenly
        SvTestHelper tester = new SvTestHelper();

        tester.logVerbose(true);
        tester.Config.RequiredAnnotations = ALL_ANNOTATIONS;

        SvVarData var1 = createTestSv("1","1","1",100,1000,-1,-1, INV,1);
        SvVarData var2 = createTestSv("2","1","1",1200,2500,1,1, INV,3);
        SvVarData var3 = createTestSv("3","1","1",2000,200,-1,1, INV,8);
        SvVarData var4 = createTestSv("4","1","1",1500,100,1,-1, DUP,12);
        SvVarData var5 = createTestSv("5","1","1",2100,2200,1,1, DEL,16);
        SvVarData var6 = createTestSv("6","1","1",500,600,-1,1, DUP,20);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.preClusteringInit();

        tester.addCopyNumberData();

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);
        assertFalse(cluster.getAnnotations().contains(CLUSTER_ANNONTATION_DM));

        tester.clearClustersAndSVs();

        // now a cluster with a set of 2 foldbacks which control the ploidies
        var1 = createTestSv("1","1","1",400,500,-1,-1, INV,1);
        var2 = createTestSv("2","1","1",5000,5200,1,1, INV,2);
        var3 = createTestSv("3","1","1",1000,1200,-1,-1, INV,2);
        var4 = createTestSv("4","1","1",2000,2200,1,1, INV,3);
        var5 = createTestSv("5","1","1",2100,2200,1,1, DEL,8);
        var6 = createTestSv("6","1","1",100,6000,-1,1, DUP,8);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.preClusteringInit();

        tester.addCopyNumberData();

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);
        assertFalse(cluster.getAnnotations().contains(CLUSTER_ANNONTATION_DM));

    }

}
