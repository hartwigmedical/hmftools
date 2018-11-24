package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createSv;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_INT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_INT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_RECIPROCAL_TRANS;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Test;

public class SyntheticDelDupTests
{
    @Test
    public void testSyntheticDelOrDupFromInvPairs()
    {
        SvTestHelper tester = new SvTestHelper();

        // create 2 INVs with varying positions to check what synthetic DEL or DUP they create

        // no overlap but 2 deletion bridges
        SvVarData var1 = createInv(tester.nextVarId(), "1", 100, 200, 1);
        SvVarData var2 = createInv(tester.nextVarId(), "1", 300, 400, -1);

        SvCluster cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_INT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), 0);

        // 1 pair of overlapping breakends

        var1 = createInv(tester.nextVarId(), "1", 100, 300, 1);
        var2 = createInv(tester.nextVarId(), "1", 200, 400, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_INT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(true));

        // test 2 DBBs but with overlapping bases less than permitted TI length
        var1 = createInv(tester.nextVarId(), "1", 100, 400, 1);
        var2 = createInv(tester.nextVarId(), "1", 90, 395, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_INT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(true));

        // test 2 overlapping breakends
        var1 = createInv(tester.nextVarId(), "1", 100, 400, 1);
        var2 = createInv(tester.nextVarId(), "1", 250, 350, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(true) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(false));

        // 3 overlapping breakends
        var1 = createInv(tester.nextVarId(), "1", 200, 400, 1);
        var2 = createInv(tester.nextVarId(), "1", 100, 350, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var1.position(true) - var2.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(false));

        // 4 overlapping breakends
        var1 = createInv(tester.nextVarId(), "1", 300, 400, 1);
        var2 = createInv(tester.nextVarId(), "1", 100, 250, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP_INT_TI);
        assertEquals(cluster.getSynDelDupLength(), var1.position(false) - var2.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(true) - var2.position(false));
    }

    @Test
    public void testSyntheticDelOrDupFromBndPairs()
    {
        SvTestHelper tester = new SvTestHelper();

        // create 2 BNDs with varying positions to check what synthetic DEL or DUP they create

        // no overlap but 2 deletion bridges - a reciprocal translation
        SvVarData var1 = createSv(tester.nextVarId(), "1", "2,", 100, 200, 1, 1, BND);
        SvVarData var2 = createSv(tester.nextVarId(), "1", "2,", 120, 220, -1, -1, BND);

        SvCluster cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markBndPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_TRANS);
        assertEquals(cluster.getSynDelDupLength(), 0);
        assertEquals(cluster.getSynDelDupTILength(), 0);

        // similar but with small overlaps
        var1 = createSv(tester.nextVarId(), "1", "2,", 100, 200, 1, 1, BND);
        var2 = createSv(tester.nextVarId(), "1", "2,", 90, 190, -1, -1, BND);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markBndPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_TRANS);
        assertEquals(cluster.getSynDelDupLength(), 0);
        assertEquals(cluster.getSynDelDupTILength(), 0);

        // single TI with a DEL
        var1 = createSv(tester.nextVarId(), "1", "2,", 100, 200, 1, 1, BND);
        var2 = createSv(tester.nextVarId(), "1", "2,", 200, 150, -1, -1, BND);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markBndPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(true) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(false));

        // 2 TIs with a DUP
        var1 = createSv(tester.nextVarId(), "1", "2,", 200, 200, 1, 1, BND);
        var2 = createSv(tester.nextVarId(), "1", "2,", 100, 150, -1, -1, BND);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markBndPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var1.position(true) - var2.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(false));
    }

    @Test
    public void testSyntheticDelOrDupFromDelsAndDups()
    {
        SvTestHelper tester = new SvTestHelper();

        // create 2 DELs next to each other
        SvVarData var1 = createDel(tester.nextVarId(), "1", 100, 200);
        SvVarData var2 = createDel(tester.nextVarId(), "1", 300, 400);

        SvCluster cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markDelDupPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_INT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var2.position(true) - var1.position(false));

        // 2 DUPs next to each other
        var1 = createDup(tester.nextVarId(), "1", 100, 200);
        var2 = createDup(tester.nextVarId(), "1", 300, 400);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markDelDupPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(true) - var1.position(false));
        assertEquals(cluster.getSynDelDupTILength(), var2.position(false) - var1.position(true));

        // 2 DUPs with an overlap
        var1 = createDup(tester.nextVarId(), "1", 100, 300);
        var2 = createDup(tester.nextVarId(), "1", 250, 400);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markDelDupPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP_INT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(true));

        // 2 DUPs with one enclosed
        var1 = createDup(tester.nextVarId(), "1", 100, 400);
        var2 = createDup(tester.nextVarId(), "1", 250, 300);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markDelDupPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(true));

        // DEL and DUP - DEL enclosed
        var1 = createDel(tester.nextVarId(), "1", 300, 350);
        var2 = createDup(tester.nextVarId(), "1", 100, 400);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markDelDupPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var1.position(true) - var2.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var2.position(false) - var1.position(false));

        // DEL and DUP - DEL overlapping DUP
        var1 = createDel(tester.nextVarId(), "1", 100, 350);
        var2 = createDup(tester.nextVarId(), "1", 200, 400);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markDelDupPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(true) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var2.position(false) - var1.position(false));

        // DEL and DUP - overlapping the other way
        var1 = createDup(tester.nextVarId(), "1", 100, 350);
        var2 = createDel(tester.nextVarId(), "1", 200, 400);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markDelDupPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(false) - var1.position(false));
        assertEquals(cluster.getSynDelDupTILength(), var2.position(true) - var1.position(true));

        // DEL before DUP no overlap
        var1 = createDel(tester.nextVarId(), "1", 100, 200);
        var2 = createDup(tester.nextVarId(), "1", 250, 400);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markDelDupPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(true) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var2.position(false) - var1.position(false));

        // DEL after DUP no overlap
        var1 = createDup(tester.nextVarId(), "1", 100, 200);
        var2 = createDel(tester.nextVarId(), "1", 250, 400);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markDelDupPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(false) - var1.position(false));
        assertEquals(cluster.getSynDelDupTILength(), var2.position(true) - var1.position(true));
    }

    private void addAndPrepareCluster(SvTestHelper tester, SvCluster cluster)
    {
        tester.clearClustersAndSVs();
        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();
        tester.Analyser.findSimpleCompleteChains();
    }

}
