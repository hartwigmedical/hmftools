package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.POLY_A_MOTIF;
import static com.hartwig.hmftools.linx.types.SvVarData.ASSEMBLY_TYPE_EQV;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSgl;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createTestSv;

import static org.junit.Assert.assertEquals;

import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

import org.junit.Test;

public class FilteringTest
{
    @Test
    public void testDuplicateBreakendFiltering()
    {
        LinxTester tester = new LinxTester();

        SvVarData var1 = createDel(tester.nextVarId(), "1", 1000, 1100);

        // equivalent breakends are kept separate
        SvVarData var2 = createSgl(tester.nextVarId(), "1", 2000, -1);

        SvVarData var3 = createSgl(tester.nextVarId(), "1", 2001, -1);

        SvVarData var4 = createSgl(tester.nextVarId(), "1", 3000, -1);

        SvVarData var5 = createSgl(tester.nextVarId(), "1", 3002, -1);
        var5.setAssemblyData(true, ASSEMBLY_TYPE_EQV);

        SvVarData var6 = createDel(tester.nextVarId(), "1", 5000, 6000);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(3, tester.getClusters().size());

        assertTrue(tester.hasClusterWithSVs(Lists.newArrayList(var1, var3, var4, var6)));

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(var2));
        assertNotNull(cluster);
        assertEquals(ResolvedType.DUP_BE, cluster.getResolvedType());

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var5));
        assertNotNull(cluster);
        assertEquals(ResolvedType.DUP_BE, cluster.getResolvedType());

        tester.clearClustersAndSVs();

        // filter out a spanning variant
        var1 = createSv(tester.nextVarId(), "1", "1", 1000, 2000, 1, -1, DEL, "1234567890123456789012345678901");

        // equivalent breakends are kept separate
        var2 = createBnd(tester.nextVarId(), "1", 999, 1, "2", 100, -1);
        var3 = createBnd(tester.nextVarId(), "1", 2001, -1, "2", 200, 1);

        var2.setAssemblyData(false, "asmn23");
        var3.setAssemblyData(false, "asmn23");

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(2, tester.getClusters().size());

        assertTrue(tester.hasClusterWithSVs(Lists.newArrayList(var2, var3)));

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var1));
        assertNotNull(cluster);
        assertEquals(ResolvedType.DUP_BE, cluster.getResolvedType());
    }

    @Test
    public void testLowVAFFiltering()
    {
        LinxTester tester = new LinxTester();

        // tester.logVerbose(true);

        SvVarData var0 = createDel(tester.nextVarId(), "1", 1000, 1100);

        // common arm rule will merge the BNDs except the one which is isolated, not poly-A and has low CNC support
        SvVarData var1 = createTestSv(tester.nextVarId(), "1", "1", 2000, 2050, -1, -1,
                INV, 1, 1, 0.1, 0.1, 0.1, "");

        SvVarData var2 = createTestSv(tester.nextVarId(), "1", "2", 10000, 10000, -1, -1,
                BND, 1, 1, 0.1, 0.1, 0.1, POLY_A_MOTIF);

        SvVarData var3 = createTestSv(tester.nextVarId(), "1", "2", 20000, 20000, -1, -1,
                BND, 1, 1, 0.1, 0.1, 0.1, "");

        SvVarData var4 = createTestSv(tester.nextVarId(), "1", "2", 30000, 38000, -1, -1,
                BND, 1, 1, 0.1, 0.1, 0.1, "");

        SvVarData var5 = createTestSv(tester.nextVarId(), "1", "2", 40000, 40000, -1, -1,
                BND, 1, 1, 0.1, 0.1, 0.1, "");

        SvVarData var6 = createTestSv(tester.nextVarId(), "1", "", 50000, -1, -1, 0,
                SGL, 1, 0, 0.1, 0, 0.1, "");

        tester.AllVariants.add(var0);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);
        tester.AllVariants.add(var6);

        tester.AllVariants.forEach(x -> x.setJcnRecalcData(0.1, 0.2));

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(5, tester.getClusters().size());

        assertTrue(tester.hasClusterWithSVs(Lists.newArrayList(var0)));
        assertTrue(tester.hasClusterWithSVs(Lists.newArrayList(var2, var4, var5)));

        SvCluster cluster = tester.findClusterWithSVs(Lists.newArrayList(var1));
        assertNotNull(cluster);
        assertEquals(ResolvedType.LOW_VAF, cluster.getResolvedType());

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var3));
        assertNotNull(cluster);
        assertEquals(ResolvedType.LOW_VAF, cluster.getResolvedType());

        cluster = tester.findClusterWithSVs(Lists.newArrayList(var6));
        assertNotNull(cluster);
        assertEquals(ResolvedType.LOW_VAF, cluster.getResolvedType());
    }

    @Test
    public void testInferredPairFiltering()
    {
        LinxTester tester = new LinxTester();

        SvVarData var0 = createTestSv(tester.nextVarId(), "1", "", 2000, 0, 1, 0,
                INF, 1, 1, 1.0, 0, 1, "");

        // common arm rule will merge the BNDs except the one which is isolated, not poly-A and has low CNC support
        SvVarData var1 = createTestSv(tester.nextVarId(), "1", "", 3000, 0, -1, 0,
                INF, 1, 1, 1.0, 0, 1, "");

        SvVarData var2 = createTestSv(tester.nextVarId(), "1", "", 4000, 0, -1, 0,
                INF, 1, 1, 1.0, 0, 1, "");

        SvVarData var3 = createTestSv(tester.nextVarId(), "1", "", 6000, 0, 1, 0,
                INF, 1, 1, 1.0, 0, 1, "");

        SvVarData var4 = createTestSv(tester.nextVarId(), "1", "", 7000, 0, -1, 0,
                INF, 1, 1, 1.0, 0, 1, "");

        SvVarData var5 = createTestSv(tester.nextVarId(), "1", "", 8000, 0, 1, 0,
                INF, 1, 1, 2.0, 0, 2, "");

        tester.AllVariants.add(var0);
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.AllVariants.add(var5);

        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(5, tester.getClusters().size());
    }

}
