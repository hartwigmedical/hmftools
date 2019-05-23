package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createBnd;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createSgl;
import static com.hartwig.hmftools.svanalysis.analysis.SvClassification.RESOLVED_TYPE_COMPLEX;
import static com.hartwig.hmftools.svanalysis.analysis.SvClassification.RESOLVED_TYPE_DEL;
import static com.hartwig.hmftools.svanalysis.analysis.SvClassification.RESOLVED_TYPE_DUP;
import static com.hartwig.hmftools.svanalysis.analysis.SvClassification.RESOLVED_TYPE_RECIPROCAL_DUP_DEL;
import static com.hartwig.hmftools.svanalysis.analysis.SvClassification.RESOLVED_TYPE_RECIPROCAL_DUP_PAIR;
import static com.hartwig.hmftools.svanalysis.analysis.SvClassification.RESOLVED_TYPE_RECIPROCAL_INV;
import static com.hartwig.hmftools.svanalysis.analysis.SvClassification.RESOLVED_TYPE_RECIPROCAL_TRANS;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_OTHER;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_TI_ONLY;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Ignore;
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

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(!cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_COMPLEX);

        // 1 pair of overlapping breakends
        var1 = createInv(tester.nextVarId(), "1", 100, 300, 1);
        var2 = createInv(tester.nextVarId(), "1", 200, 400, -1);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(false) - var2.position(true));

        assertEquals(1, cluster.getArmClusters().size());
        assertEquals(ARM_CL_DSB, cluster.getArmClusters().get(0).getType());
        assertEquals(1, cluster.getArmClusters().get(0).getTICount());

        // test 2 DSBs but with an overlapping end less than permitted TI length
        var1 = createInv(tester.nextVarId(), "1", 100, 400, 1);
        var2 = createInv(tester.nextVarId(), "1", 90, 1000, -1);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(false) - var2.position(true));


        // test 2 overlapping breakends
        var1 = createInv(tester.nextVarId(), "1", 100, 10400, 1);
        var2 = createInv(tester.nextVarId(), "1", 250, 10350, -1);
        var1.setAssemblyData(false, "asmb12");
        var2.setAssemblyData(false, "asmb12");

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var2.position(true) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(false) - var2.position(false));

        assertEquals(2, cluster.getArmClusters().size());
        assertEquals(ARM_CL_DSB, cluster.getArmClusters().get(0).getType());
        assertEquals(ARM_CL_TI_ONLY, cluster.getArmClusters().get(1).getType());

        // test 2 overlapping breakends but where a pair of breakend form an overlapping DB
        var1 = createInv(tester.nextVarId(), "1", 500, 1000, 1);
        var2 = createInv(tester.nextVarId(), "1", 480, 2000, -1);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(false) - var2.position(true));

        // 3 overlapping breakends
        var1 = createInv(tester.nextVarId(), "1", 200, 10400, 1);
        var2 = createInv(tester.nextVarId(), "1", 100, 10350, -1);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP);
        assertEquals(cluster.getSyntheticLength(), var1.position(true) - var2.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(false) - var2.position(false));

        assertEquals(2, cluster.getArmClusters().size());
        assertEquals(ARM_CL_COMPLEX_OTHER, cluster.getArmClusters().get(0).getType());
        assertEquals(ARM_CL_TI_ONLY, cluster.getArmClusters().get(1).getType());

        // 4 overlapping breakends - also forming 2 facing foldbacks
        var1 = createInv(tester.nextVarId(), "1", 300, 400, 1);
        var2 = createInv(tester.nextVarId(), "1", 100, 250, -1);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP);
        assertEquals(cluster.getSyntheticLength(), var1.position(false) - var2.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(true) - var2.position(false));

        assertEquals(1, cluster.getArmClusters().size());
        assertEquals(ARM_CL_COMPLEX_FOLDBACK, cluster.getArmClusters().get(0).getType());
    }

    @Test
    public void testSyntheticDelDupFromBndPairs()
    {
        SvTestHelper tester = new SvTestHelper();

        // create 2 BNDs with varying positions to check what synthetic DEL or DUP they create

        // single TI with a DEL
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 100, 1, "2", 200, 1);
        SvVarData var2 = createBnd(tester.nextVarId(), "1", 200, -1, "2", 150, -1);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var2.position(true) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(false) - var2.position(false));

        // 2 TIs with a DUP
        var1 = createBnd(tester.nextVarId(), "1", 200, 1, "2", 200, 1);
        var2 = createBnd(tester.nextVarId(), "1", 100, -1, "2", 150, -1);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP);
        assertEquals(cluster.getSyntheticLength(), var1.position(true) - var2.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(false) - var2.position(false));
    }

    @Test
    public void testSyntheticDelDupFromDelsAndDups()
    {
        SvTestHelper tester = new SvTestHelper();

        // 2 DUPs next to each other
        SvVarData var1 = createDup(tester.nextVarId(), "1", 100, 200);
        SvVarData var2 = createDup(tester.nextVarId(), "1", 300, 400);

        var1.getTIAssemblies(true).add("asmb1");
        var2.getTIAssemblies(false).add("asmb1");

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var2.position(true) - var1.position(false));
        assertEquals(cluster.getSyntheticTILength(), var2.position(false) - var1.position(true));

        // 2 DUPs with an overlap
        var1 = createDup(tester.nextVarId(), "1", 100, 300);
        var2 = createDup(tester.nextVarId(), "1", 250, 400);

        var1.getTIAssemblies(false).add("asmb1");
        var2.getTIAssemblies(true).add("asmb1");

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP);
        assertEquals(cluster.getSyntheticLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(false) - var2.position(true));

        // 2 DUPs with one enclosed
        var1 = createDup(tester.nextVarId(), "1", 100, 400);
        var2 = createDup(tester.nextVarId(), "1", 250, 300);

        var1.getTIAssemblies(false).add("asmb1");
        var2.getTIAssemblies(true).add("asmb1");

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP);
        assertEquals(cluster.getSyntheticLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(false) - var2.position(true));

        // DEL and DUP - DEL enclosed
        var1 = createDel(tester.nextVarId(), "1", 300, 350);
        var2 = createDup(tester.nextVarId(), "1", 100, 400);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP);
        assertEquals(cluster.getSyntheticLength(), var1.position(true) - var2.position(true));
        assertEquals(cluster.getSyntheticTILength(), var2.position(false) - var1.position(false));

        // DEL and DUP - DEL overlapping DUP
        var1 = createDel(tester.nextVarId(), "1", 100, 350);
        var2 = createDup(tester.nextVarId(), "1", 200, 400);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var2.position(true) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var2.position(false) - var1.position(false));

        // DEL and DUP - overlapping the other way
        var1 = createDup(tester.nextVarId(), "1", 100, 350);
        var2 = createDel(tester.nextVarId(), "1", 200, 400);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var2.position(false) - var1.position(false));
        assertEquals(cluster.getSyntheticTILength(), var2.position(true) - var1.position(true));

        // DEL before DUP no overlap
        var1 = createDel(tester.nextVarId(), "1", 100, 200);
        var2 = createDup(tester.nextVarId(), "1", 250, 400);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var2.position(true) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var2.position(false) - var1.position(false));

        // DEL after DUP no overlap
        var1 = createDup(tester.nextVarId(), "1", 100, 200);
        var2 = createDel(tester.nextVarId(), "1", 250, 400);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var2.position(false) - var1.position(false));
        assertEquals(cluster.getSyntheticTILength(), var2.position(true) - var1.position(true));
    }

    @Test
    public void testSyntheticChainedDelsAndDups()
    {
        SvTestHelper tester = new SvTestHelper();

        SvVarData var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 1000, -1);
        SvVarData var2 = createBnd(tester.nextVarId(), "2", 1200, 1, "3", 2000, 1);
        SvVarData var3 = createBnd(tester.nextVarId(), "1", 5000, -1, "3", 1600, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DEL);
        assertEquals(cluster.getSyntheticLength(), var3.position(true) - var1.position(true));
        assertEquals(400, cluster.getSyntheticTILength()); // takes the longest

        tester.clearClustersAndSVs();

        // add one more SV to extend the chain into a DUP - the INV comes in too close for a TI with the first SV
        SvVarData var4 = createDup(tester.nextVarId(), "1", 980, 5600);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_DUP);
        assertTrue(cluster.isSyntheticType());
        assertEquals(cluster.getSyntheticLength(), var1.position(true) - var4.position(true));
        assertEquals(600, cluster.getSyntheticTILength()); // takes average

        tester.clearClustersAndSVs();

        // special case of the reciprocal DUP which is more likely to not be chained through a long TI, and so is split
        var1 = createBnd(tester.nextVarId(), "1", 20000, 1, "3", 1000, -1);
        var2 = createBnd(tester.nextVarId(), "2", 1000, -1, "3", 1500, 1);

        var3 = createBnd(tester.nextVarId(), "1", 1000, -1, "3", 2000, -1);
        var4 = createBnd(tester.nextVarId(), "2", 20000, 1, "3", 2500, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);
        assertEquals(2, cluster.getChains().size());

        assertTrue(cluster.isResolved());
        assertEquals(RESOLVED_TYPE_RECIPROCAL_DUP_PAIR, cluster.getResolvedType());
        assertTrue(!cluster.isSyntheticType());
        assertEquals(0, cluster.getSyntheticLength());
        assertEquals(500, cluster.getSyntheticTILength());

    }

    @Test
    public void testSyntheticReciprocalInversions()
    {
        // test reciprocal inversion defined as a TI greater than 50% of the synthetic length and internal to a DEL
        SvTestHelper tester = new SvTestHelper();

        SvVarData var1 = createInv(tester.nextVarId(), "1", 100, 1000, 1);
        SvVarData var2 = createInv(tester.nextVarId(), "1", 90, 990, -1);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_INV);
        assertEquals(cluster.getSyntheticLength(), var2.position(false) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var1.position(false) - var2.position(true));

        // again but with a longer chain - initially with the DB too long
        var1 = createBnd(tester.nextVarId(), "1", 1000, 1, "2", 1500, 1);
        var2 = createBnd(tester.nextVarId(), "1", 20000, -1, "2", 1000, -1);
        SvVarData var3 = createDel(tester.nextVarId(), "1", 22000, 25000);

        tester.clearClustersAndSVs();
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_DUP_DEL);

        var2 = createBnd(tester.nextVarId(), "1", 2000, -1, "2", 1000, -1);
        var3 = createBnd(tester.nextVarId(), "1", 21000, 1, "2", 5000, -1);
        SvVarData var4 = createBnd(tester.nextVarId(), "1", 21500, -1, "2", 5500, 1);

        tester.clearClustersAndSVs();
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_INV);
        assertEquals(cluster.getSyntheticLength(), var4.position(true) - var1.position(true));
        assertEquals(cluster.getSyntheticTILength(), var3.position(true) - var2.position(true));

    }

    @Test
    public void testSyntheticReciprocalTranslocations()
    {
        SvTestHelper tester = new SvTestHelper();

        // 2 BNDs - no overlap but 2 deletion bridges - a reciprocal translation
        SvVarData var1 = createBnd(tester.nextVarId(), "1", 100, 1, "2", 200, 1);
        SvVarData var2 = createBnd(tester.nextVarId(), "1", 120, -1, "2", 220, -1);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_TRANS);
        assertEquals(cluster.getSyntheticLength(), 0);
        assertEquals(cluster.getSyntheticTILength(), 0);

        // similar but with small overlaps
        var1 = createBnd(tester.nextVarId(), "1", 100, 1, "2", 200, 1);
        var2 = createBnd(tester.nextVarId(), "1", 90, -1, "2", 190, -1);

        addAndCluster(tester, var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.Analyser.getClusters().get(0);

        assertTrue(cluster.isResolved());
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_TRANS);
        assertEquals(cluster.getSyntheticLength(), 0);
        assertEquals(cluster.getSyntheticTILength(), 0);

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
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_TRANS);
        assertEquals(200, cluster.getSyntheticTILength()); // takes average


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
        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_RECIPROCAL_TRANS);
        assertEquals(600, cluster.getSyntheticTILength()); // takes average
    }

    @Ignore
    @Test
    public void testSyntheticDelOrDupExcludedCases()
    {
        SvTestHelper tester = new SvTestHelper();

        // non-TI breakends cannot be on different arms - makes them inconsistent any way

        // DUP enclosing DEL
        SvVarData var1 = createDup(tester.nextVarId(), "1", 100, 150000000);
        SvVarData var2 = createDel(tester.nextVarId(), "1", 130000000, 140000000);

        addAndCluster(tester, var1, var2);

        assertEquals(2, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        // assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_NONE);

        // 2 INVs but with a complex breakend in the middle of the TI
        var1 = createInv(tester.nextVarId(), "1", 100, 50000, 1);
        var2 = createInv(tester.nextVarId(), "1", 300, 20000, -1);

        SvVarData sgl = createSgl(tester.nextVarId(), "1", 35000, 1, false);

        tester.clearClustersAndSVs();
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(sgl);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(2, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.getResolvedType() == RESOLVED_TYPE_COMPLEX);
    }

    private void addAndCluster(SvTestHelper tester, SvVarData var1, SvVarData var2)
    {
        tester.clearClustersAndSVs();
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();
    }

}
