package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createBnd;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createSgl;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createSv;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_OTHER;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_ISOLATED_BE;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_TI_ONLY;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.getArmClusterData;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_COMPLEX;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_INT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_INT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_NONE;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_RECIPROCAL_INV;
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
    public void testSyntheticDelDupFromInvPairs()
    {
        SvTestHelper tester = new SvTestHelper();

        // create 2 INVs with varying positions to check what synthetic DEL or DUP they create

        // no overlap but 2 deletion bridges - doesn't make anything
        SvVarData var1 = createInv(tester.nextVarId(), "1", 100, 200, 1);
        SvVarData var2 = createInv(tester.nextVarId(), "1", 300, 400, -1);

        SvCluster cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(!cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_COMPLEX);

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

        cluster.buildArmClusters();
        assertEquals(1, cluster.getArmClusters().size());
        assertEquals(ARM_CL_DSB, cluster.getArmClusters().get(0).getType());
        assertEquals(1, cluster.getArmClusters().get(0).getTICount());

        // test 2 DSBs but with an overlapping end less than permitted TI length
        var1 = createInv(tester.nextVarId(), "1", 100, 400, 1);
        var2 = createInv(tester.nextVarId(), "1", 90, 1000, -1);

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
        var1 = createInv(tester.nextVarId(), "1", 100, 10400, 1);
        var2 = createInv(tester.nextVarId(), "1", 250, 10350, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(true) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(false));

        cluster.buildArmClusters();
        assertEquals(2, cluster.getArmClusters().size());
        assertEquals(ARM_CL_DSB, cluster.getArmClusters().get(0).getType());
        assertEquals(ARM_CL_TI_ONLY, cluster.getArmClusters().get(1).getType());

        // test 2 overlapping breakends but where a pair of breakend form an overlapping DB
        var1 = createInv(tester.nextVarId(), "1", 500, 1000, 1);
        var2 = createInv(tester.nextVarId(), "1", 480, 2000, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_INT_TI);
        assertEquals(cluster.getSynDelDupLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(true));

        // 3 overlapping breakends
        var1 = createInv(tester.nextVarId(), "1", 200, 10400, 1);
        var2 = createInv(tester.nextVarId(), "1", 100, 10350, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP_EXT_TI);
        assertEquals(cluster.getSynDelDupLength(), var1.position(true) - var2.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(false));

        cluster.buildArmClusters();
        assertEquals(2, cluster.getArmClusters().size());
        assertEquals(ARM_CL_COMPLEX_OTHER, cluster.getArmClusters().get(0).getType());
        assertEquals(ARM_CL_TI_ONLY, cluster.getArmClusters().get(1).getType());

        // 4 overlapping breakends - also forming 2 facing foldbacks
        var1 = createInv(tester.nextVarId(), "1", 300, 400, 1);
        var2 = createInv(tester.nextVarId(), "1", 100, 250, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);
        tester.Analyser.markFoldbacks();

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP_INT_TI);
        assertEquals(cluster.getSynDelDupLength(), var1.position(false) - var2.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(true) - var2.position(false));

        cluster.buildArmClusters();
        assertEquals(1, cluster.getArmClusters().size());
        assertEquals(ARM_CL_COMPLEX_FOLDBACK, cluster.getArmClusters().get(0).getType());


        // test reciprocal inversion defined as a TI greater than 50% of the synthetic length and internal to a DEL
        var1 = createInv(tester.nextVarId(), "1", 100, 1000, 1);
        var2 = createInv(tester.nextVarId(), "1", 90, 990, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_INV);
        assertEquals(cluster.getSynDelDupLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSynDelDupTILength(), var1.position(false) - var2.position(true));

    }

    @Test
    public void testSyntheticDelDupFromBndPairs()
    {
        SvTestHelper tester = new SvTestHelper();

        // create 2 BNDs with varying positions to check what synthetic DEL or DUP they create

        // no overlap but 2 deletion bridges - a reciprocal translation
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 100, 1, "2", 200, 1);
        SvVarData var2 = createBnd(tester.nextVarId(), "1", 120, -1, "2", 220, -1);

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
        var1 = createBnd(tester.nextVarId(), "1", 100, 1, "2", 200, 1);
        var2 = createBnd(tester.nextVarId(), "1", 90, -1, "2", 190, -1);

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
        var1 = createBnd(tester.nextVarId(), "1", 100, 1, "2", 200, 1);
        var2 = createBnd(tester.nextVarId(), "1", 200, -1, "2", 150, -1);

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
        var1 = createBnd(tester.nextVarId(), "1", 200, 1, "2", 200, 1);
        var2 = createBnd(tester.nextVarId(), "1", 100, -1, "2", 150, -1);

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
    public void testSyntheticDelDupFromDelsAndDups()
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

        var1.getTempInsertionAssemblies(true).add("asmb1");
        var2.getTempInsertionAssemblies(false).add("asmb1");

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

        var1.getTempInsertionAssemblies(false).add("asmb1");
        var2.getTempInsertionAssemblies(true).add("asmb1");

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

        var1.getTempInsertionAssemblies(false).add("asmb1");
        var2.getTempInsertionAssemblies(true).add("asmb1");

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

    @Test
    public void testSyntheticDelOrDupExcludedCases()
    {
        SvTestHelper tester = new SvTestHelper();

        // non-TI breakends cannot be on different arms - makes them inconsistent any way

        // DUP enclosing DEL
        SvVarData var1 = createDup(tester.nextVarId(), "1", 100, 150000000);
        SvVarData var2 = createDel(tester.nextVarId(), "1", 130000000, 140000000);

        SvCluster cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markDelDupPairTypes(cluster);

        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_NONE);

        // 2 INVs but with a complex breakend in the middle of the TI
        var1 = createInv(tester.nextVarId(), "1", 100, 50000, 1);
        var2 = createInv(tester.nextVarId(), "1", 300, 20000, -1);

        cluster = new SvCluster(0);
        cluster.addVariant(var1);
        cluster.addVariant(var2);

        addAndPrepareCluster(tester, cluster);
        tester.ClusteringMethods.markInversionPairTypes(cluster);
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL_EXT_TI);

        tester.clearClustersAndSVs();
        cluster.setResolved(false, RESOLVED_TYPE_NONE);

        SvCluster cluster2 = new SvCluster(1);
        SvVarData sgl = createSgl(tester.nextVarId(), "1", 35000, 1, false);
        cluster2.addVariant(sgl);

        tester.addClusterAndSVs(cluster);
        tester.addClusterAndSVs(cluster2);
        tester.preClusteringInit();
        tester.Analyser.findLimitedChains();

        tester.ClusteringMethods.markInversionPairTypes(cluster);

        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_COMPLEX);
    }

    private void addAndPrepareCluster(SvTestHelper tester, SvCluster cluster)
    {
        tester.clearClustersAndSVs();
        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();
        tester.Analyser.findLimitedChains();
    }

}
