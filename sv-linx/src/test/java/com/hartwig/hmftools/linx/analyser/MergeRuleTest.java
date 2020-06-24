package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.CONSEC_BREAKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.FOLDBACKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.HOM_LOSS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.LOH_CHAIN;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.LONG_DEL_DUP_INV;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.MAJOR_ALLELE_JCN;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.OVERLAP_FOLDBACKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringReason.TI_JCN_MATCH;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createIns;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSgl;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createTestSv;
import static com.hartwig.hmftools.linx.types.SvVarData.ASSEMBLY_TYPE_EQV;

import static org.junit.Assert.assertEquals;

import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import com.hartwig.hmftools.linx.utils.LinxTester;

public class MergeRuleTest
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
        // simple clustered SVs are split out in the final routine except for the SVs joined by assembly
        LinxTester tester = new LinxTester();

        tester.clearClustersAndSVs();

        SvVarData del1 = createDel(tester.nextVarId(), "2", 1000, 3000);
        tester.AllVariants.add(del1);

        SvVarData del2 = createDel(tester.nextVarId(), "2", 2000, 4000);
        tester.AllVariants.add(del2);

        SvVarData del3 = createDel(tester.nextVarId(), "2", 3500, 8000);
        tester.AllVariants.add(del3);

        // another group remains a cluster due to assembly linking the variants
        SvVarData var1 = createDel(tester.nextVarId(), "1", 20000, 25000);
        SvVarData var2 = createDel(tester.nextVarId(), "1", 26500, 50000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        var1.setAssemblyData(false, "asmb12");
        var2.setAssemblyData(true, "asmb12");

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(4, tester.getClusters().size());

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(var1, var2));
        assertTrue(cluster != null);
        assertEquals(ResolvedType.SIMPLE_GRP, cluster.getResolvedType());
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
    public void testFoldbackMerge()
    {
        LinxTester tester = new LinxTester();

        // 2 clusters with foldbacks on the same arm are merged
        SvVarData inv1 = createInv(tester.nextVarId(), "1", 100, 200, -1);
        tester.AllVariants.add(inv1);

        SvVarData inv2 = createInv(tester.nextVarId(), "1", 20000, 20100, 1);
        tester.AllVariants.add(inv2);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.getClusters().size());
        assertTrue(inv1.hasClusterReason(FOLDBACKS));
        assertTrue(inv2.hasClusterReason(FOLDBACKS));
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
        assertTrue(cluster != null);
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
        assertTrue(cluster != null);
        assertTrue(cluster.hasClusterReason(HOM_LOSS));
        assertTrue(cluster.hasClusterReason(MAJOR_ALLELE_JCN));

        /*
        cluster = tester.findClusterWithSVs(Lists.newArrayList(var2, var3));
        assertTrue(cluster != null);
        assertTrue(cluster.hasClusterReason(CR_HOM_LOSS));

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var4, var5));
        assertTrue(cluster != null);
        assertTrue(cluster.hasClusterReason(CR_HOM_LOSS));
        */


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
        assertTrue(cluster != null);
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
        assertTrue(cluster != null);
        assertTrue(cluster.hasClusterReason(HOM_LOSS));
    }

    @Test
    public void testLohResolvingClusterMerge()
    {
        // merge clusters based on their SVs being required to stop a SV chaining through an LOH event
        LinxTester tester = new LinxTester();

        SvVarData var1 = createDel(1, "1", 10000, 60000);
        SvVarData var2 = createDup(2, "1", 20000, 50000);
        SvVarData var3 = createBnd(3, "1", 30000, 1, "2", 10000, 1);
        SvVarData var4 = createBnd(4, "1", 40000, -1, "2", 20000, 1);
        SvVarData var5 = createDel(5, "1", 70000, 120000);
        SvVarData var6 = createDup(6, "1", 80000, 110000);
        SvVarData var7 = createBnd(7, "1", 90000, 1, "2", 80000, 1);
        SvVarData var8 = createBnd(8, "1", 100000, -1, "2", 90000, -1);
        SvVarData var9 = createDup(9, "2", 30000, 60000);

        // does run unto the DUP in the LOH but isn't clustered since is simple
        SvVarData var10 = createDel(10, "2", 40000, 50000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);
        tester.AllVariants.add(var8);
        tester.AllVariants.add(var9);
        tester.AllVariants.add(var10);

        List<LohEvent> lohData = tester.CnDataLoader.getLohData();

        lohData.add(new LohEvent("1", 10000, 20000, "DEL", "DUP", 1, var1.id(), var2.id()));

        lohData.add(new LohEvent("1", 50000, 60000,"DUP", "DEL", 1, var2.id(), var1.id()));

        lohData.add(new LohEvent("1", 110000, 120000, "DUP", "DEL", 1, var6.id(), var5.id()));

        lohData.add(new LohEvent("2", 20000, 30000,"BND", "DUP", 1, var4.id(), var9.id()));

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(2, tester.Analyser.getClusters().size());
        assertTrue(var10.getCluster().getSvCount() == 1);

        assertTrue(var2.hasClusterReason(LOH_CHAIN));
        assertTrue(var2.getClusterReason().contains(String.valueOf(var4.id())));
        assertTrue(var4.hasClusterReason(LOH_CHAIN));
        assertTrue(var4.getClusterReason().contains(String.valueOf(var2.id())));
    }

    @Test
    public void testConsistentBreakendOverlapMerge()
    {
        LinxTester tester = new LinxTester();

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

        // eqv breakend will be ignored
        SvVarData var4 = createSgl(tester.nextVarId(), "1", 80000, -1);
        var4.setAssemblyData(true, ASSEMBLY_TYPE_EQV);
        allVariants.add(var4);

        // now some SV which will be overlapped by the consecutive breakends
        SvVarData overlap1 = createBnd(tester.nextVarId(), "1", 20000, -1, "4", 200, 1);
        allVariants.add(overlap1);

        SvVarData overlap2 = createSgl(tester.nextVarId(), "1", 40000, -1);
        allVariants.add(overlap2);

        tester.AllVariants.addAll(allVariants);
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(4, tester.getClusters().size());

        SvCluster mainCluster = tester.findClusterWithSVs(Lists.newArrayList(consec1, consec2, consec3, overlap1, overlap2));

        if(mainCluster == null)
            assertTrue(false);

        assertTrue(overlap1.hasClusterReason(CONSEC_BREAKS));
        assertTrue(overlap2.hasClusterReason(CONSEC_BREAKS));
    }

    @Test
    public void testFoldbacksStraddlingBreakendMerge()
    {
        LinxTester tester = new LinxTester();

        // a cluster with 2 foldbacks and an unclustered breakend in between

        // foldbacks clustered due to the foldback rule
        SvVarData var1 = createInv(tester.nextVarId(), "1", 1000, 2000, -1);
        SvVarData var2 = createInv(tester.nextVarId(), "1", 100000, 101000, -1);

        // simple del not clustered
        SvVarData var3 = createDel(tester.nextVarId(), "1", 20000, 21000);

        SvVarData var4 = createBnd(tester.nextVarId(), "1", 40000, -1, "3", 100, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.isFoldback());
        assertTrue(var2.isFoldback());

        assertEquals(2, tester.getClusters().size());

        SvCluster mainCluster = tester.findClusterWithSVs(Lists.newArrayList(var1, var2, var4));

        assertTrue(mainCluster != null);

        assertTrue(var4.hasClusterReason(OVERLAP_FOLDBACKS));
        assertTrue(mainCluster.hasClusterReason(OVERLAP_FOLDBACKS));

        // don't merge if both foldbacks face away
        tester.clearClustersAndSVs();

        var1 = createInv(tester.nextVarId(), "1", 1000, 2000, 1);

        // var2 as before

        var3 = createBnd(tester.nextVarId(), "1", 40000, -1, "3", 100, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.isFoldback());
        assertTrue(var2.isFoldback());

        assertEquals(2, tester.getClusters().size());
    }

    @Test
    public void testDistancePloidyLinkMatchMerge()
    {
        LinxTester tester = new LinxTester();

        // set a high major allele ploidy on the other chromatid to stop the MAP merge rule kicking in
        tester.setNonClusterAllelePloidies(4, 1);

        // 2 distance clusters which can form a distant ploidy-matching merge

        // foldbacks clustered due to the foldback rule
        SvVarData var0 = createTestSv(tester.nextVarId(), "1", "0", 500, -1, -1, -1, SGL, 2);
        SvVarData var1 = createTestSv(tester.nextVarId(), "1", "2", 1000, 100, 1, -1, BND, 2);
        SvVarData var2 = createTestSv(tester.nextVarId(), "1", "2", 2000, 2000, -1, 1, BND, 2);

        // simple del not clustered
        SvVarData var3 = createDel(tester.nextVarId(), "1", 20000, 21000);

        // pair of BNDs resolved to a simple type
        SvVarData var4 = createTestSv(tester.nextVarId(), "1", "3", 50000, 100, 1, -1, BND, 1);
        SvVarData var5 = createTestSv(tester.nextVarId(), "1", "3", 51000, 200, -1, 1, BND, 1);

        SvVarData var6 = createTestSv(tester.nextVarId(), "1", "1", 100000, 110000, 1, 1, INV, 2);

        tester.AllVariants.add(var0);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(3, tester.getClusters().size());

        assertTrue(tester.hasClusterWithSVs(Lists.newArrayList(var4, var5)));

        SvCluster mainCluster = tester.findClusterWithSVs(Lists.newArrayList(var0, var1, var2, var6));

        assertTrue(mainCluster != null);

        assertTrue(var6.hasClusterReason(TI_JCN_MATCH));
        assertTrue(mainCluster.hasClusterReason(TI_JCN_MATCH));

    }

}
