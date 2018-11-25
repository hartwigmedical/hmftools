package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createBnd;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDup;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createIns;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createSgl;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.ASSEMBLY_TYPE_EQV;

import static org.junit.Assert.assertEquals;

import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Test;

public class MergeRuleTests
{

    @Test
    public void testProximityMerge()
    {
        SvTestHelper tester = new SvTestHelper();

        // basic proximity clustering
        SvVarData var1 = createDel(tester.nextVarId(), "1", 1000, 1100);
        tester.AllVariants.add(var1);

        SvVarData var2 = createDup(tester.nextVarId(), "1", 2000, 2100);
        tester.AllVariants.add(var2);

        SvVarData var3 = createIns(tester.nextVarId(), "1", 3000, 3100);
        tester.AllVariants.add(var3);

        SvVarData var4 = createInv(tester.nextVarId(), "1", 4000, 4100, 1);
        tester.AllVariants.add(var4);

        SvVarData var5 = createSgl(tester.nextVarId(), "1", 5000, -1, false);
        tester.AllVariants.add(var5);

        SvVarData var6 = createBnd(tester.nextVarId(), "1", 6000, -1, "2", 1000, 1);
        tester.AllVariants.add(var6);

        // equivalent breakends are kept separate
        SvVarData var7 = createSgl(tester.nextVarId(), "1", 7000, -1, false);
        tester.AllVariants.add(var7);

        SvVarData var8 = createSgl(tester.nextVarId(), "1", 7000, -1, false);
        tester.AllVariants.add(var8);

        SvVarData var9 = createSgl(tester.nextVarId(), "1", 7100, -1, false);
        var9.setAssemblyData(true, ASSEMBLY_TYPE_EQV);
        tester.AllVariants.add(var9);

        tester.preClusteringInit();

        mergeOnProximity(tester);

        assertEquals(tester.getClusters().size(), 4);
        assertEquals(tester.getClusters().get(0).getCount(), 6);
        assertTrue(tester.getClusters().get(1).getSVs().contains(var7));
        assertTrue(tester.getClusters().get(2).getSVs().contains(var8));
        assertTrue(tester.getClusters().get(3).getSVs().contains(var9));

        // non-overlapping DELs can merge, whereas overlapping DELs are split out
        tester.clearClustersAndSVs();

        SvVarData del1 = createDel(tester.nextVarId(), "2", 1000, 3000);
        tester.AllVariants.add(del1);

        SvVarData del2 = createDel(tester.nextVarId(), "2", 2000, 4000);
        tester.AllVariants.add(del2);

        SvVarData del3 = createDel(tester.nextVarId(), "2", 3500, 8000);
        tester.AllVariants.add(del3);

        SvVarData del4 = createDel(tester.nextVarId(), "2", 6000, 10000);
        tester.AllVariants.add(del4);

        mergeOnProximity(tester);

        assertEquals(tester.getClusters().size(), 2);
        assertTrue(tester.getClusters().get(0).getSVs().contains(del1));
        assertTrue(tester.getClusters().get(0).getSVs().contains(del3));
        assertTrue(tester.getClusters().get(0).getSVs().contains(del4));
        assertTrue(tester.getClusters().get(1).getSVs().contains(del2));

        // test exclusion of low-quality SVs
        tester.clearClustersAndSVs();
    }

    private void mergeOnProximity(SvTestHelper tester)
    {
        tester.ClusteringMethods.clusterByBaseDistance(tester.AllVariants, tester.Analyser.getClusters());
    }
}
