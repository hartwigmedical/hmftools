package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createBnd;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createSgl;
import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createTestSv;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.ALL_ANNOTATIONS;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_FOLDBACK;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_FOLDBACK_DSB;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_ISOLATED_BE;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_TI_ONLY;
import static com.hartwig.hmftools.linx.types.SvArmCluster.getArmClusterData;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNONTATION_DM;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Ignore;
import org.junit.Test;

public class AnnotationTest
{
    @Test
    public void testLocalTopology()
    {
        SvTestHelper tester = new SvTestHelper();
        // tester.logVerbose(true);

        // create a set of SVs which exhibit the various types of local topology
        final SvVarData var1 = createInv("1", "1", 101000, 107000, 1);

        // remote TI from BND, with other ends forming a DB and an isolated BE
        final SvVarData var2 = createBnd("2", "1", 149000, 1, "2", 1000, -1);
        final SvVarData var3 = createBnd("3", "1", 250000, 1, "2", 1100, 1);

        // DUP forms a DB with the next foldback and a DB with the previous BND
        final SvVarData var4 = createDup("4", "1", 150000, 200500);
        final SvVarData var5 = createInv("5", "1", 201000, 204000, -1);

        // complex foldback
        final SvVarData var6 = createInv("6", "1", 301000, 304000, 1);
        final SvVarData var7 = createDel("7", "1", 305000, 306000);
        final SvVarData var8 = createInv("8", "1", 310000, 315000, -1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);
        tester.AllVariants.add(var8);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.isFoldback());
        assertTrue(var5.isFoldback());
        assertTrue(var6.isFoldback());
        assertTrue(var8.isFoldback());

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(0);

        final int[] armClusterData = getArmClusterData(cluster);

        assertEquals(6, cluster.getArmClusters().size());
        assertEquals(1, armClusterData[ARM_CL_ISOLATED_BE]);
        assertEquals(1, armClusterData[ARM_CL_TI_ONLY]);
        assertEquals(1, armClusterData[ARM_CL_DSB]);
        assertEquals(1, armClusterData[ARM_CL_FOLDBACK]);
        assertEquals(1, armClusterData[ARM_CL_FOLDBACK_DSB]);
        assertEquals(1, armClusterData[ARM_CL_COMPLEX_FOLDBACK]);
    }

    @Test
    public void testUnderClusteredFoldbackAnnotations()
    {
        SvTestHelper tester = new SvTestHelper();
        // tester.logVerbose(true);
        tester.Config.RequiredAnnotations = ALL_ANNOTATIONS;

        double[] chrCopyNumbers = {2.0, 2.0, 2.0};
        tester.ClusteringMethods.getChrCopyNumberMap().put("1", chrCopyNumbers);
        tester.ClusteringMethods.getChrCopyNumberMap().put("2", chrCopyNumbers);

        final SvVarData var1 = createInv("0", "1", 101000, 104000, -1);

        // straddling BNDs
        final SvVarData var2 = createBnd("1", "1", 1000, 1, "2", 1000, -1);
        final SvVarData var3 = createBnd("2", "1", 200000, 1, "2", 1100, 1);

        // single other cluster - will be merged in because it faces the foldback
        final SvVarData var4 = createSgl("3", "1", 150000, 1, false);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertTrue(var1.isFoldback());
        assertEquals(2, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.Analyser.getClusters().get(1);
        assertTrue(cluster.getSVs().contains(var4));

        tester.clearClustersAndSVs();

        final SvVarData var5 = createInv(tester.nextVarId(), "1", 1000, 2000, -1);
        final SvVarData var6 = createInv(tester.nextVarId(), "1", 3000, 5000, -1);
        final SvVarData var7 = createInv(tester.nextVarId(), "1", 4000, 6000, 1);
        final SvVarData var8 = createInv(tester.nextVarId(), "1", 7000, 8000, 1);

        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);
        tester.AllVariants.add(var7);
        tester.AllVariants.add(var8);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
    }

}
