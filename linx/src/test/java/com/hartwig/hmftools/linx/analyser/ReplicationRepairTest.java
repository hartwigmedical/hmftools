package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_REP_REPAIR;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class ReplicationRepairTest
{
    @Test
    public void testReplicationRepairOneBreakDupInversions()
    {
        // a single section if replicated and connected back in incorrect orientations, looking like INV - FB INV - DEL
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createInv(  1, "1", 1000, 10000, 1);
        final SvVarData var2 = createDup(  2, "1", 1100, 9500);
        final SvVarData var3 = createInv(3, "1", 1500, 10100, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        tester.Analyser.annotateClusters();

        assertEquals(1, cluster.getChains().size());
        assertTrue(cluster.hasAnnotation(CLUSTER_ANNOT_REP_REPAIR));
    }

    @Test
    public void testReplicationRepairOneBreakFoldbackInversion()
    {
        // a single section if replicated and connected back in incorrect orientations, looking like INV - FB INV - DEL
        LinxTester tester = new LinxTester();

        final SvVarData var1 = createInv(  1, "1", 1000, 10000, 1);
        final SvVarData var2 = createInv(2, "1", 1100, 1500, -1);
        final SvVarData var3 = createDel(  3, "1", 9500, 10100);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertTrue(var2.isFoldback());

        assertEquals(1, tester.Analyser.getClusters().size());
        final SvCluster cluster = tester.Analyser.getClusters().get(0);

        tester.Analyser.annotateClusters();

        assertEquals(1, cluster.getChains().size());
        assertTrue(cluster.hasAnnotation(CLUSTER_ANNOT_REP_REPAIR));
    }
}
