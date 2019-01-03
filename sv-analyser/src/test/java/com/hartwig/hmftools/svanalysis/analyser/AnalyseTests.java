package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createSgl;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createTestSv;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Test;

public class AnalyseTests
{
    @Test
    public void testArmClusterAnalysis()
    {
        SvTestHelper tester = new SvTestHelper();

        // form a cluster with multiple DSBs and subsequent replication and breakage


        // 2 x DSBs
        // main section: F F B F B B F F B B
        // 0 -> 1000 -> 1100 -> 1200 -> 1400 , 2 - 1 - 2 - 3 - 4
        // 1400 -> 5400 -> 5600 -> 8600, 4 - 3 - 4 - 3
        // 8600 -> 9000 -> 13020 -> 14200, 3 - 2 - 3 - 4
        // 14200 -> 16200 -> 18200 -> 18400 -> 19000, 4 - 3 - 2 - 1 - 2

        // SvVarData var3 = createDel(tester.nextVarId(), "1", "1", 1000, 1100, 1, -1, DEL, 2, 2, 1, 1, 1, "");
        tester.AllVariants.add(createTestSv(tester.nextVarId(), "1", "1", 1000, 1100, 1, -1, DEL, 2, 2, 1, 1, 1, ""));

        // SvVarData sgl = createSgl(tester.nextVarId(), "1", 1000, -1, false);
        // tester.AllVariants.add(createSgl(tester.nextVarId(), "1", 1200, -1, false));
        // tester.AllVariants.add(createInv(tester.nextVarId(), "1", 1200, 13020, -1, 3, 1));
        tester.AllVariants.add(createTestSv(tester.nextVarId(), "1", "1", 1200, 13020, -1, -1, INV, 3, 3, 1, 1, 1, ""));

        // tester.AllVariants.add(createSgl(tester.nextVarId(), "1", 1400, -1, false));
        // tester.AllVariants.add(createDup(tester.nextVarId(), "1", 1400, 9000));
        tester.AllVariants.add(createTestSv(tester.nextVarId(), "1", "1", 1400, 9000, -1, 1, DUP, 4, 3, 1, 1, 1, ""));

        // tester.AllVariants.add(createSgl(tester.nextVarId(), "1", 5400, 1, false));
        tester.AllVariants.add(createTestSv(tester.nextVarId(), "1", "1", 5400, -1, 1, -1, SGL, 4, 0, 1, 0, 1, ""));
        // tester.AllVariants.add(createSgl(tester.nextVarId(), "1", 5600, -1, false));
        tester.AllVariants.add(createTestSv(tester.nextVarId(), "1", "1", 5600, -1, -1, -1, SGL, 4, 0, 1, 0, 1, ""));
        // tester.AllVariants.add(createSgl(tester.nextVarId(), "1", 8600, 1, false));
        tester.AllVariants.add(createTestSv(tester.nextVarId(), "1", "1", 8600, -1, 1, -1, SGL, 4, 0, 1, 0, 1, ""));

        //tester.AllVariants.add(createSgl(tester.nextVarId(), "1", 9000, 1, false));
        // tester.AllVariants.add(createSgl(tester.nextVarId(), "1", 13020, -1, false));

        // tester.AllVariants.add(createSgl(tester.nextVarId(), "1", 14200, -1, false));
        // tester.AllVariants.add(createSgl(tester.nextVarId(), "1", 16200, 1, false));
        tester.AllVariants.add(createTestSv(tester.nextVarId(), "1", "1", 16200, -1, 1, -1, SGL, 4, 0, 1, 0, 1, ""));
        // tester.AllVariants.add(createDup(tester.nextVarId(), "1", 14200, 18200));
        tester.AllVariants.add(createTestSv(tester.nextVarId(), "1", "1", 14200, 18200, -1, 1, DUP, 4, 3, 1, 1, 1, ""));
        // tester.AllVariants.add(createSgl(tester.nextVarId(), "1", 18200, 1, false));
        // tester.AllVariants.add(createDel(tester.nextVarId(), "1", 18400, 19000));
        tester.AllVariants.add(createTestSv(tester.nextVarId(), "1", "1", 18400, 19000, 1, -1, DEL, 2, 2, 1, 1, 1, ""));

        tester.preClusteringInit();

        tester.mergeOnProximity();
        assertEquals(tester.getClusters().size(), 1);

        tester.Analyser.clusterAndAnalyse();




    }
}