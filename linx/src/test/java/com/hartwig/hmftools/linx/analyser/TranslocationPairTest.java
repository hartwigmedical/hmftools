package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.linx.types.ResolvedType.DEL_TI;
import static com.hartwig.hmftools.linx.types.ResolvedType.DUP_TI;
import static com.hartwig.hmftools.linx.types.ResolvedType.PAIR_OTHER;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS_DEL_DUP;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS_DUPS;
import static com.hartwig.hmftools.linx.types.ResolvedType.UNBAL_TRANS_TI;
import static com.hartwig.hmftools.linx.types.LinxConstants.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;

import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import com.hartwig.hmftools.linx.utils.LinxTester;

public class TranslocationPairTest
{

    @Test
    public void testReciprocalTranslocations()
    {
        LinxTester tester = new LinxTester();

        // 2 BNDs - no overlap but 2 deletion bridges - a reciprocal translation
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 100, 1, "2", 200, 1);
        SvVarData var2 = createBnd(tester.nextVarId(), "1", 120, -1, "2", 220, -1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RECIP_TRANS);

        // similar but with small overlaps
        var1 = createBnd(tester.nextVarId(), "1", 100, 1, "2", 200, 1);
        var2 = createBnd(tester.nextVarId(), "1", 90, -1, "2", 190, -1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RECIP_TRANS);

        tester.clearClustersAndSVs();

        // from a chain and a BND
        var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 1000, -1);
        var2 = createBnd(tester.nextVarId(), "2", 1200, 1, "3", 2000, -1);
        SvVarData var3 = createBnd(tester.nextVarId(), "1", 4000, -1, "3", 1600, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RECIP_TRANS);

        tester.clearClustersAndSVs();

        // 2 chains
        var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 1000, -1);
        var2 = createBnd(tester.nextVarId(), "2", 1200, 1, "3", 2000, -1);

        var3 = createBnd(tester.nextVarId(), "1", 4000, -1, "4", 10000, 1);
        SvVarData var4 = createBnd(tester.nextVarId(), "3", 1600, 1, "4", 9400, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RECIP_TRANS);
        assertEquals(2, cluster.getChains().size());
        assertFalse(cluster.getChains().stream().anyMatch(x -> x.getLinkedPairs().stream().anyMatch(y -> y.baseLength() > SHORT_TI_LENGTH)));
    }

    @Test
    public void testFacingDelTranslocations()
    {
        LinxTester tester = new LinxTester();

        SvVarData var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 50001, -1);
        SvVarData var2 = createBnd(tester.nextVarId(), "1", 5000, -1, "2", 150000, 1);

        tester.CnDataLoader.getLohData().add(new LohEvent(var1.chromosome(false), 1, var1.position(false),
                SegmentSupport.TELOMERE.toString(), var1.typeStr(), 1, LohEvent.CN_DATA_NO_SV, var1.id()));

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DEL_TI);

        // and now as a RECIP_DEL_DUP
        var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 50001, -1);
        var2 = createBnd(tester.nextVarId(), "1", 5000, -1, "2", 150000, 1);

        tester.CnDataLoader.getLohData().clear();

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RECIP_TRANS_DEL_DUP);
        assertTrue(cluster.getChains().isEmpty());
    }

    @Test
    public void testFacingDupTranslocations()
    {
        LinxTester tester = new LinxTester();

        SvVarData var1 = createBnd(tester.nextVarId(), "1", 1000, -1, "2", 50000, -1);
        SvVarData var2 = createBnd(tester.nextVarId(), "1", 200000, 1, "2", 150000, 1);

        tester.CnDataLoader.getLohData().add(new LohEvent(var1.chromosome(false), 1, var1.position(false),
                SegmentSupport.TELOMERE.toString(), var1.typeStr(), 1, LohEvent.CN_DATA_NO_SV, var1.id()));

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(!cluster.isResolved());
        assertTrue(cluster.getResolvedType() == DUP_TI);

        // and now as a RECIP_DUPS
        var1 = createBnd(tester.nextVarId(), "1", 1000, -1, "2", 50000, -1);
        var2 = createBnd(tester.nextVarId(), "1", 200000, 1, "2", 150000, 1);

        tester.CnDataLoader.getLohData().clear();

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertTrue(!cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RECIP_TRANS_DUPS);
        assertTrue(cluster.getChains().isEmpty());

        // chained version
        var1 = createBnd(tester.nextVarId(), "1", 1000, -1, "3", 100, -1);
        var2 = createBnd(tester.nextVarId(), "2", 50000, -1, "3", 200, 1);

        SvVarData var3 = createBnd(tester.nextVarId(), "1", 200000, 1, "4", 100, -1);
        SvVarData var4 = createBnd(tester.nextVarId(), "2", 150000, 1, "4", 200, 1);

        tester.clearClustersAndSVs();
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        cluster = tester.getClusters().get(0);
        assertTrue(!cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RECIP_TRANS_DUPS);
        assertEquals(2, cluster.getChains().size());
    }

    @Test
    public void testPairOtherTranslocations()
    {
        LinxTester tester = new LinxTester();

        // 2 BNDs with 3 arms
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 1000, -1, "2", 1000, -1);
        SvVarData var2 = createBnd(tester.nextVarId(), "2", 10000, 1, "3", 1000, 1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(!cluster.isResolved());
        assertTrue(cluster.getResolvedType() == UNBAL_TRANS_TI);

        // BNDs not forming any TIs and not a balanced translation aren't classified
        tester.clearClustersAndSVs();

        var1 = createBnd(tester.nextVarId(), "1", 1000, -1, "2", 1000, -1);
        var2 = createBnd(tester.nextVarId(), "2", 5000, -1, "3", 1000, 1);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertTrue(!cluster.isResolved());
        assertTrue(cluster.getResolvedType() == PAIR_OTHER);

    }

}
