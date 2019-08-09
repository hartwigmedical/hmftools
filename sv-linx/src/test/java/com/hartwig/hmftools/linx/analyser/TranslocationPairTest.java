package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.linx.analysis.SvClassification.getSyntheticLength;
import static com.hartwig.hmftools.linx.analysis.SvClassification.getSyntheticTiLength;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import static Utils.SvTestRoutines.createBnd;

import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import Utils.LinxTester;

public class TranslocationPairTest
{

    @Test
    public void testSyntheticReciprocalTranslocations()
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
    }

}
