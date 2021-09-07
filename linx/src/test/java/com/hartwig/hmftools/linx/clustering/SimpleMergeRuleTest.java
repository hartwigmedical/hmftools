package com.hartwig.hmftools.linx.clustering;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.HIGH_JCN;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.HOM_LOSS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.LONG_DEL_DUP_INV;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.MAJOR_ALLELE_JCN;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createIns;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSgl;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createTestSv;

import static org.junit.Assert.assertEquals;

import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class SimpleMergeRuleTest
{
    @Test
    public void testProximityMerge()
    {
        LinxTester tester = new LinxTester();

        // basic proximity clustering with 2 duplicate breakend SGLs excluded
        SvVarData var1 = createDel(tester.nextVarId(), "1", 1000, 1100);
        tester.AllVariants.add(var1);

        SvVarData var2 = createDup(tester.nextVarId(), "1", 2000, 2100);
        tester.AllVariants.add(var2);

        SvVarData var3 = createIns(tester.nextVarId(), "1", 3000, 3100);
        tester.AllVariants.add(var3);

        SvVarData var4 = createInv(tester.nextVarId(), "1", 4000, 4100, 1);
        tester.AllVariants.add(var4);

        SvVarData var5 = createSgl(tester.nextVarId(), "1", 5000, -1);
        tester.AllVariants.add(var5);

        SvVarData var6 = createBnd(tester.nextVarId(), "1", 6000, -1, "5", 1000, 1);
        tester.AllVariants.add(var6);

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
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.getClusters().size());

        assertTrue(tester.hasClusterWithSVs(Lists.newArrayList(var1, var2, var3, var4, var5, var6, var10,
                var11, var12, var13)));
    }

    @Test
    public void testSimpleSVsDemerge()
    {
        // simple clustered SVs of 2 SVs are split out in the final routine, 3+ are kept as COMPLEX
        LinxTester tester = new LinxTester();

        tester.clearClustersAndSVs();

        SvVarData del1 = createDel(tester.nextVarId(), "2", 1000, 3000);
        tester.AllVariants.add(del1);

        SvVarData del2 = createDel(tester.nextVarId(), "2", 2000, 4000);
        tester.AllVariants.add(del2);

        SvVarData del3 = createDel(tester.nextVarId(), "2", 3500, 8000);
        tester.AllVariants.add(del3);

        SvVarData var1 = createDel(tester.nextVarId(), "1", 20000, 25000);
        SvVarData var2 = createDel(tester.nextVarId(), "1", 26500, 50000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        var1.setAssemblyData(false, "asmb12");
        var2.setAssemblyData(true, "asmb12");

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(3, tester.getClusters().size());

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(del1, del2, del3));
        assertNotNull(cluster);
        assertEquals(ResolvedType.COMPLEX, cluster.getResolvedType());
    }

    @Test
    public void testLongDelDupInvMerge()
    {
        LinxTester tester = new LinxTester();

        // INV, DUP and DEL overlapping or enclosed are merged if within long merge distance and can form a TI
        SvVarData var1 = createInv(tester.nextVarId(), "1", 1000, 200000, -1);
        SvVarData var2 = createInv(tester.nextVarId(), "1", 20000, 220000, -1);

        // not merged since face the same way
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();
        assertEquals(2, tester.getClusters().size());

        tester.clearClustersAndSVs();

        // now an INV which faces the others will cause them both to be clustered in
        SvVarData var3 = createInv(tester.nextVarId(), "1", 40000, 240000, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();
        assertEquals(1, tester.getClusters().size());

        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(var1.hasClusterReason(LONG_DEL_DUP_INV));
        assertTrue(cluster.hasClusterReason(LONG_DEL_DUP_INV));

        // test again with some very distance but overlapping SVs ignored, and others in a DB included
        tester.clearClustersAndSVs();

        // now an INV which faces the others will cause them both to be clustered in
        var1 = createDel(tester.nextVarId(), "1", 10000, 30000000);

        // INV faces same way as DEL on start so not clustered and other breakends are too far apart
        var2 = createInv(tester.nextVarId(), "1", 500000, 10000000, 1);

        // DUP faces but is too far away
        var3 = createDup(tester.nextVarId(), "1", 20000000, 40000000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();
        assertEquals(3, tester.getClusters().size());

        tester.clearClustersAndSVs();

        // a group of 3 DELs even if far apart will be merged if they overlap with with another group of DUPs 3 or more times
        var1 = createDup(tester.nextVarId(), "1", 10000000, 30001000);
        var2 = createDup(tester.nextVarId(), "1", 30000000, 50001000);
        var3 = createDup(tester.nextVarId(), "1", 50000000, 70001000);
        SvVarData var4 = createInv(tester.nextVarId(), "1", 70000000, 71000000, -1);

        SvVarData var5 = createDel(tester.nextVarId(), "1", 1000000, 20000000);
        SvVarData var6 = createDel(tester.nextVarId(), "1", 20001000, 40000000);
        SvVarData var7 = createDel(tester.nextVarId(), "1", 40001000, 60000000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();
        assertEquals(1, tester.getClusters().size());

        cluster = tester.getClusters().get(0);

        assertTrue(tester.AllVariants.stream().anyMatch(x -> x.hasClusterReason(LONG_DEL_DUP_INV)));
        assertTrue(cluster.hasClusterReason(LONG_DEL_DUP_INV));
    }

    @Test
    public void testHighFacingJcnMerge()
    {
        LinxTester tester = new LinxTester();

        // test 2 overlapping DUPs, too far for proximity, are merged
        SvVarData var1 = createDel(1, "1", 1000, 2000);

        SvVarData var2 = createTestSv(2, "1", "1", 10000,40000, -1, 1, DUP,  6);
        SvVarData var3 = createTestSv(3, "1", "1", 30000,70000, -1, 1, DUP,  6);

        // to prevent a single group being dissolved
        SvVarData var4 = createTestSv(4, "1", "1", 72000,80000, -1, -1, INV,  2);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4));
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(2, tester.Analyser.getClusters().size());

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(var2, var3, var4));
        assertNotNull(cluster);

        assertTrue(cluster.hasClusterReason(HIGH_JCN));

        tester.clearClustersAndSVs();

        // a more complicated version

        var1 = createDel(1, "1", 1000, 2000);

        var2 = createTestSv(2, "1", "", 10000,0, -1, 0, SGL,  40);

        // first SGL is too small a CN drop to satisfy the rule
        var3 = createTestSv(3, "1", "", 30000,0, 1, 0, SGL,  6);

        // simple SV2 ignored
        var4 = createDel(4, "1", 40000, 41000);
        SvVarData var5 = createDup(5, "1", 50000, 51000);

        // next SGL raises the CN again
        SvVarData var6 = createTestSv(6, "1", "", 60000,0, -1, 0, SGL,  10);

        // both the next 2 satisfy the drop in CN criteria
        SvVarData var7 = createTestSv(7, "1", "", 70000,0, 1, 0, SGL,  36);
        SvVarData var8 = createTestSv(8, "1", "", 80000,0, 1, 0, SGL,  8);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4, var5, var6, var7, var8));
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(4, tester.Analyser.getClusters().size());

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var2, var3, var6, var7, var8));
        assertNotNull(cluster);

        assertTrue(cluster.hasClusterReason(HIGH_JCN));
        assertTrue(cluster.hasClusterReason(MAJOR_ALLELE_JCN));

        tester.clearClustersAndSVs();

        // now 2 clusters each formed from this rule but kept separated from each other
        var1 = createDel(1, "1", 1000, 2000);

        var2 = createTestSv(2, "1", "", 10000,0, -1, 0, SGL,  40);
        var3 = createTestSv(3, "1", "", 30000,0, 1, 0, SGL,  30);
        var4 = createTestSv(4, "1", "", 40000,0, 1, 0, SGL,  8);

        var5 = createDup(5, "1", 50000, 51000);
        var6 = createDel(6, "1", 60000, 61000);

        var7 = createTestSv(7, "1", "1", 70000,100000, -1, 1, DUP,  10);
        var8 = createTestSv(8, "1", "1", 90000,120000, -1, 1, DUP,  10);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4, var5, var6, var7, var8));
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(5, tester.Analyser.getClusters().size());

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var2, var3, var4));
        assertNotNull(cluster);
        assertTrue(cluster.hasClusterReason(HIGH_JCN));

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var7, var8));
        assertNotNull(cluster);
        assertTrue(cluster.hasClusterReason(HIGH_JCN));
    }

    @Test
    public void testLohHomLossEventMerge()
    {
        LinxTester tester = new LinxTester();

        // scenario 1: LOH containing all clustered HOM-loss events should also be clustered
        SvVarData var1 = createBnd(1, "1", 1000, 1, "2", 100, 1);
        SvVarData var2 = createBnd(2, "1", 100000, -1, "3", 100, 1);

        // 2x hom-loss events both clustered
        SvVarData var3 = createDel(3, "1", 6500, 6600);
        SvVarData var4 = createBnd(4, "1", 20000, 1, "5", 200, 1);
        SvVarData var5 = createBnd(5, "1", 22000, -1, "5", 100, -1);

        LohEvent lohEvent = new LohEvent("1", 1000, 100000, "BND", "BND", 1, var1.id(), var2.id());

        tester.CnDataLoader.getLohData().add(lohEvent);

        List<HomLossEvent> homLossData = tester.CnDataLoader.getHomLossData();

        homLossData.add(new HomLossEvent(var3.chromosome(true), var3.position(true), var3.position(false),
                var3.typeStr(), var3.typeStr(), var3.id(), var3.id()));

        homLossData.add(new HomLossEvent(var4.chromosome(true), var4.position(true), var5.position(true),
                var4.typeStr(), var5.typeStr(), var4.id(), var5.id()));

        lohEvent.addHomLossEvents(homLossData);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(3, tester.Analyser.getClusters().size());

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(var1, var2));
        assertNotNull(cluster);
        assertTrue(cluster.hasClusterReason(HOM_LOSS));

        // scenario 2: multiple hom-loss events clustered because LOH is clustered
        tester.clearClustersAndSVs();

        var1 = createDel(1, "1", 1000, 100000);
        var2 = createBnd(2, "1", 10000, 1, "2", 100, 1);
        var3 = createBnd(3, "1", 20000, -1, "3", 100, 1);
        var4 = createBnd(4, "1", 30000, 1, "4", 100, 1);
        var5 = createBnd(5, "1", 40000, -1, "5", 100, 1);

        tester.CnDataLoader.getLohData().clear();

        lohEvent = new LohEvent(var1.chromosome(true), var1.position(true), var1.position(false),
                "DEL", "DEL", 1, var1.id(), var1.id());

        tester.CnDataLoader.getLohData().add(lohEvent);

        homLossData.clear();

        homLossData.add(new HomLossEvent(var2.chromosome(true), var2.position(true), var3.position(true),
                var2.typeStr(), var3.typeStr(), var2.id(), var3.id()));

        homLossData.add(new HomLossEvent(var4.chromosome(true), var4.position(true), var5.position(true),
                var4.typeStr(), var5.typeStr(), var4.id(), var5.id()));

        lohEvent.addHomLossEvents(homLossData);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        // the 2 hom-loss deletes are clustered with each other due to the major allele ploidy rule
        assertEquals(2, tester.Analyser.getClusters().size());

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var2, var3, var4, var5));
        assertNotNull(cluster);
        assertTrue(cluster.hasClusterReason(HOM_LOSS));
        assertTrue(cluster.hasClusterReason(MAJOR_ALLELE_JCN));

        // scenario 3: hom-loss event overlaps a LOH
        tester.clearClustersAndSVs();

        var1 = createDel(1, "1", 10000, 30000);
        var2 = createBnd(2, "1", 20000, 1, "2", 100, 1);
        var3 = createBnd(3, "1", 40000, -1, "3", 100, 1);

        tester.CnDataLoader.getLohData().clear();

        lohEvent = new LohEvent(var1.chromosome(true), var1.position(true), var3.position(true),
                var1.typeStr(), var3.typeStr(), 1,var1.id(), var3.id());

        tester.CnDataLoader.getLohData().add(lohEvent);

        homLossData.clear();

        homLossData.add(new HomLossEvent(var2.chromosome(true), var2.position(true), var1.position(false),
                var2.typeStr(), var1.typeStr(), var2.id(), var1.id()));

        lohEvent.addHomLossEvents(homLossData);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var2, var3));
        assertNotNull(cluster);
        assertTrue(cluster.hasClusterReason(HOM_LOSS));

        // again but with more SVs involved
        tester.clearClustersAndSVs();

        var1 = createDel(1, "1", 10000, 30000);
        var2 = createDel(2, "1", 40000, 60000);
        var3 = createBnd(3, "1", 20000, 1, "2", 100, 1);
        var4 = createBnd(4, "1", 50000, -1, "3", 100, 1);

        tester.CnDataLoader.getLohData().clear();

        lohEvent = new LohEvent(var1.chromosome(true), var1.position(true), var2.position(false),
                var1.typeStr(), var2.typeStr(), 1, var1.id(), var2.id());

        tester.CnDataLoader.getLohData().add(lohEvent);

        homLossData.clear();

        homLossData.add(new HomLossEvent(var2.chromosome(true), var3.position(true), var1.position(false),
                var3.typeStr(), var1.typeStr(), var3.id(), var1.id()));

        homLossData.add(new HomLossEvent(var2.chromosome(true), var2.position(true), var4.position(true),
                var2.typeStr(), var4.typeStr(), var2.id(), var4.id()));

        lohEvent.addHomLossEvents(homLossData);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var3, var4));
        assertNotNull(cluster);
        assertTrue(cluster.hasClusterReason(HOM_LOSS));
    }


}
