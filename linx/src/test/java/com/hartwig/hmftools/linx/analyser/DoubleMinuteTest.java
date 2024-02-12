package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.analysis.DoubleMinuteData.SEG_DATA_COUNT;
import static com.hartwig.hmftools.linx.types.ResolvedType.COMPLEX;
import static com.hartwig.hmftools.linx.types.ResolvedType.DOUBLE_MINUTE;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS_DEL_DUP;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createTestSv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.analysis.DoubleMinuteData;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Ignore;
import org.junit.Test;

public class DoubleMinuteTest
{
    @Test
    public void testSimpleDupDM()
    {
        LinxTester tester = new LinxTester();

        // need to put another SV before it to set the background copy number for this chromatid
        final SvVarData var1 = createTestSv(1,"1","1",500,600,1,-1, DEL,1);

        // first a simple DUP
        final SvVarData dup = createTestSv(2,"1","1",50000,55000,-1,1, DUP,10);
        dup.setJcnRecalcData(8, 12);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(dup);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(2, tester.Analyser.getClusters().size());

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(dup));
        assertTrue(cluster != null);
        assertEquals(1, cluster.getSvCount());
        assertTrue(cluster.getSVs().contains(dup));
        assertEquals(DOUBLE_MINUTE, cluster.getResolvedType());

        assertTrue(cluster.getChains().size() == 1);
        assertEquals(1, cluster.getChains().get(0).getLinkCount());
        assertTrue(cluster.getChains().get(0).isClosedLoop());
    }

    @Test
    public void testChainedDM()
    {
        // form a DM from 2 INVs
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createTestSv(1,"1","1",1000,6000,-1,-1, INV,8);
        final SvVarData var2 = createTestSv(2,"1","1",3000,8000,1,1, INV,8);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(var1, var2));
        assertTrue(cluster != null);
        assertEquals(DOUBLE_MINUTE, cluster.getResolvedType());
        assertTrue(cluster.getDoubleMinuteSVs().contains(var1));
        assertTrue(cluster.getDoubleMinuteSVs().contains(var2));

        assertEquals(1, cluster.getChains().size());
        SvChain chain = cluster.getChains().get(0);
        assertEquals(2, chain.getLinkCount());
        assertEquals(2, chain.getSvCount());
        assertTrue(chain.isClosedLoop());
    }

    @Test
    public void testInvalidOpenBreakendDMs()
    {
        // a cluster is not a DM if it has too high a ratio of open breakends, with simple DELs excluded
        LinxTester tester = new LinxTester();

        // the main DM structure, with open breakends
        SvVarData var1 = createTestSv(tester.nextVarId(),"1","2",1001,100,-1,1, BND,8);
        SvVarData var2 = createTestSv(tester.nextVarId(),"1","2",20000,200,1,-1, BND,8);
        SvVarData var3 = createTestSv(tester.nextVarId(),"1","1",2000,10001,1,-1, DEL,8);
        SvVarData var4 = createTestSv(tester.nextVarId(),"1","1",11000,16000,1,-1, DEL,8);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4));
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getDoubleMinuteFinder().getDoubleMinutes().size());
        DoubleMinuteData dmData = tester.Analyser.getDoubleMinuteFinder().getDoubleMinutes().get(0);

        assertEquals(2, dmData.SimpleDels);
        assertEquals(2, dmData.OpenBreakends);
        assertEquals(2, dmData.ClosedBreakends);

        assertEquals(1, tester.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);
        assertEquals(RECIP_TRANS_DEL_DUP, cluster.getResolvedType());
    }

    @Test
    public void testInternalExternalBreakendDMs()
    {
        // a cluster is not a DM if it has too many breakends going from within a segment to outside, but  short TIs are excluded
        LinxTester tester = new LinxTester();

        // tester.logVerbose(true);

        // the main DM structure
        SvVarData var1 = createTestSv(tester.nextVarId(),"1","2",1000,1000,-1,-1, BND,8);
        SvVarData var2 = createTestSv(tester.nextVarId(),"1","2",10000,10000,1,1, BND,8);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(var1, var2));
        assertTrue(cluster != null);
        assertEquals(DOUBLE_MINUTE, cluster.getResolvedType());
        assertTrue(cluster.getDoubleMinuteSVs().contains(var1));
        assertTrue(cluster.getDoubleMinuteSVs().contains(var2));

        tester.clearClustersAndSVs();

        var1 = createTestSv(tester.nextVarId(),"1","2",1000,1000,-1,-1, BND,8);
        var2 = createTestSv(tester.nextVarId(),"1","2",10000,10000,1,1, BND,8);

        // add short TIs and test again - invalidating the DM
        SvVarData var3 = createTestSv(tester.nextVarId(),"1","3",2000,100,1,-1, BND,1.5);
        SvVarData var4 = createTestSv(tester.nextVarId(),"1","3",2100,200,-1,1, BND,1.5);

        SvVarData var5 = createTestSv(tester.nextVarId(),"1","4",3000,100,-1,1, BND,1.5);
        SvVarData var6 = createTestSv(tester.nextVarId(),"1","4",3100,200,1,-1, BND,1.5);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4, var5, var6));
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getDoubleMinuteFinder().getDoubleMinutes().size());
        DoubleMinuteData dmData = tester.Analyser.getDoubleMinuteFinder().getDoubleMinutes().get(0);

        assertEquals(1, dmData.ChainIntExtData.size());
        assertEquals(4, (int)dmData.ChainIntExtData.get(dmData.Chains.get(0).id())[SEG_DATA_COUNT]);

        assertEquals(1, tester.getClusters().size());
        cluster = tester.getClusters().get(0);
        assertTrue(cluster != null);
        assertEquals(COMPLEX, cluster.getResolvedType());

        tester.clearClustersAndSVs();

        var1 = createTestSv(tester.nextVarId(),"1","2",1000,1000,-1,-1, BND,8);
        var2 = createTestSv(tester.nextVarId(),"1","2",10000,10000,1,1, BND,8);

        // this time assembled links will be ignored
        var3 = createTestSv(tester.nextVarId(),"1","3",2000,100,1,-1, BND,1.5);
        var4 = createTestSv(tester.nextVarId(),"1","3",2100,200,-1,1, BND,1.5);
        var3.setAssemblyData(false, "asmb34");
        var4.setAssemblyData(false, "asmb34");

        var5 = createTestSv(tester.nextVarId(),"1","4",3000,100,-1,1, BND,1.5);
        var6 = createTestSv(tester.nextVarId(),"1","4",3100,200,1,-1, BND,1.5);
        var5.setAssemblyData(true, "asmb56");
        var5.setAssemblyData(true, "asmb56");

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4, var5, var6));
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getDoubleMinuteFinder().getDoubleMinutes().size());
        dmData = tester.Analyser.getDoubleMinuteFinder().getDoubleMinutes().get(0);

        assertTrue(dmData.ChainIntExtData.isEmpty());

        assertEquals(1, tester.getClusters().size());
        cluster = tester.getClusters().get(0);
        assertTrue(cluster != null);
        assertEquals(COMPLEX, cluster.getResolvedType());
        assertTrue(cluster.getAnnotations().contains(CLUSTER_ANNOT_DM));
        assertEquals(2, cluster.getDoubleMinuteSVs().size());
    }

    @Test
    public void testInvalidChainedFoldbackDMs()
    {
        // if a foldback splits a chain, it should be the foldback's JCN which is tested
        LinxTester tester = new LinxTester();

        // tester.logVerbose(true);

        // the main DM structure
        SvVarData var1 = createTestSv(tester.nextVarId(),"1","1",1000,1100,-1,-1, INV,8);
        SvVarData var2 = createTestSv(tester.nextVarId(),"1","2",5000,1000,1,-1, BND,16);
        SvVarData var3 = createTestSv(tester.nextVarId(),"2","2",5100,5000,1,1, INV,8);

        // an SV which will invalidate the DM
        SvVarData var4 = createTestSv(tester.nextVarId(),"1","3",2000,1000,1,1, BND,2);
        SvVarData var5 = createTestSv(tester.nextVarId(),"1","4",3000,1000,-1,1, BND,2);
        SvVarData var6 = createTestSv(tester.nextVarId(),"1","5",4000,1000,1,1, BND,2);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4, var5, var6));
        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);
        assertTrue(cluster != null);
        assertEquals(COMPLEX, cluster.getResolvedType());
    }

    @Test
    public void testMultipleDMsInCluster()
    {
        LinxTester tester = new LinxTester();

        SvVarData var1 = createTestSv(1,"1","1",1000,10000,-1,1, DUP,8);
        SvVarData var2 = createTestSv(2,"1","1",2000,3000,1,-1, DEL,1);
        var1.setAssemblyData(true, "asmb12");
        var2.setAssemblyData(true, "asmb12");

        SvVarData var3 = createTestSv(3,"1","1",12000,12500,-1,-1, INV,8);
        SvVarData var4 = createTestSv(4,"1","1",20000,20500,1,1, INV,8);
        SvVarData var5 = createTestSv(5,"1","1",13000,14000,1,-1, DEL,1);
        SvVarData var6 = createTestSv(6,"1","2",15000,100,1,-1, BND,1);
        SvVarData var7 = createTestSv(7,"1","2",16000,200,-1,1, BND,1);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4, var5, var6, var7));

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster != null);

        assertEquals(COMPLEX, cluster.getResolvedType());
        assertTrue(cluster.getAnnotations().contains(CLUSTER_ANNOT_DM));

        assertTrue(cluster.getDoubleMinuteSVs().contains(var1));
        assertTrue(cluster.getDoubleMinuteSVs().contains(var3));
        assertTrue(cluster.getDoubleMinuteSVs().contains(var4));

        assertEquals(2, cluster.getChains().stream().filter(x -> x.isDoubleMinute()).count());

        // test again but with 1 chain valid, and 2 not
        tester.clearClustersAndSVs();

        var1 = createTestSv(tester.nextVarId(),"1","1",25000,35000,-1,1, DUP,8);

        var2 = createTestSv(tester.nextVarId(),"1","2",100,1100,-1,1, BND,8);
        var3 = createTestSv(tester.nextVarId(),"1","2",10000,1000,1,-1, BND,8);
        var4 = createTestSv(tester.nextVarId(),"1","1",1000,20000,1,1, INV,1.5);
        var5 = createTestSv(tester.nextVarId(),"1","1",2000,21000,-1,-1, INV,1.5);
        var6 = createTestSv(tester.nextVarId(),"1","1",3000,22000,1,1, INV,1.5);

        // this DUP is too short
        var7 = createTestSv(tester.nextVarId(),"1","1",36000,36500,-1,1, DUP,8);

        tester.AllVariants.addAll(Lists.newArrayList(var1, var2, var3, var4, var5, var6, var7));

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);
        assertTrue(cluster != null);

        assertEquals(COMPLEX, cluster.getResolvedType());
        assertTrue(cluster.getAnnotations().contains(CLUSTER_ANNOT_DM));

        assertEquals(1, cluster.getDoubleMinuteSVs().size());
        assertTrue(cluster.getDoubleMinuteSVs().contains(var1));

        assertEquals(1, cluster.getChains().stream().filter(x -> x.isDoubleMinute()).count());
    }

    @Ignore
    @Test
    public void testChainedDMWithMinors()
    {
        // NOTE: set to ignored since cannot form into a single chain

        // form a DM from 3 chained SVs, with some other SVs in the cluster having a lower JCN
        LinxTester tester = new LinxTester();

        // 1 s10100 -> 6 e10600-10500s -> 4 s11500-10100e -> 3 s10200-12000e -> 5 s12100-12200e -> 2 e12500-11200s -> 1 e11000

        final SvVarData var1 = createTestSv(1,"1","1",10100,11000,-1,-1, INV,8);
        final SvVarData var2 = createTestSv(2,"1","1",11200,12500,1,1, INV,8);
        final SvVarData var3 = createTestSv(3,"1","2",12000,10200,-1,1, BND,8);
        final SvVarData var4 = createTestSv(4,"1","2",11500,10100,1,-1, BND,8);
        final SvVarData var5 = createTestSv(5,"1","1",12100,12200,1,-1, DEL,2);
        final SvVarData var6 = createTestSv(6,"1","1",10500,10600,-1,1, DUP,2);

        // unrelated SVs at either end of the cluster
        final SvVarData other1 = createTestSv(7,"1","1",100,200,1,-1, DEL,1);
        final SvVarData other2 = createTestSv(8,"1","1",20000,20100,1,-1, DEL,1);
        final SvVarData other3 = createTestSv(9,"2","2",20100,20100,1,-1, DEL,1);

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

        tester.Analyser.clusterAndAnalyse();

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(var1, var2, var3, var4, var5, var6));
        assertTrue(cluster != null);
        assertTrue(cluster.hasAnnotation(CLUSTER_ANNOT_DM));
        assertTrue(cluster.getDoubleMinuteSVs().contains(var1));
        assertTrue(cluster.getDoubleMinuteSVs().contains(var2));
        assertTrue(cluster.getDoubleMinuteSVs().contains(var3));
        assertTrue(cluster.getDoubleMinuteSVs().contains(var4));

        assertTrue(cluster.getChains().size() == 1);
        SvChain chain = cluster.getChains().get(0);
        assertEquals(18, chain.getLinkCount());
        assertEquals(6, chain.getSvCount());
    }

}
