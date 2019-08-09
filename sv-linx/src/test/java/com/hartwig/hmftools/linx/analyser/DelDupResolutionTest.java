package com.hartwig.hmftools.linx.analyser;


import static com.hartwig.hmftools.linx.analysis.SvClassification.getSyntheticGapLength;
import static com.hartwig.hmftools.linx.analysis.SvClassification.getSyntheticLength;
import static com.hartwig.hmftools.linx.analysis.SvClassification.getSyntheticTiLength;
import static com.hartwig.hmftools.linx.types.ResolvedType.DEL_TI;
import static com.hartwig.hmftools.linx.types.ResolvedType.DUP_TI;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import static com.hartwig.hmftools.linx.utils.SvTestRoutines.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestRoutines.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestRoutines.createDup;

import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import com.hartwig.hmftools.linx.utils.LinxTester;

public class DelDupResolutionTest
{

    @Test
    public void testEnclosedDelInDup()
    {
        LinxTester tester = new LinxTester();

        SvVarData var1 = createDup(tester.nextVarId(), "1", 1000, 50000);
        SvVarData var2 = createDel(tester.nextVarId(), "1", 5000, 40000);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(cluster.getResolvedType() == DUP_TI);
        assertEquals(var1.position(false) - var2.position(false), getSyntheticLength(cluster));
        assertEquals(var2.position(true) - var1.position(true), getSyntheticTiLength(cluster));
        assertEquals(var2.position(false) - var2.position(true), getSyntheticGapLength(cluster));
    }

    @Test
    public void testOverlappingDelDup()
    {
        LinxTester tester = new LinxTester();

        SvVarData var1 = createDup(tester.nextVarId(), "1", 1000, 40000);
        SvVarData var2 = createDel(tester.nextVarId(), "1", 36000, 60000);

        tester.addAndCluster(var1, var2);

        assertEquals(1, tester.Analyser.getClusters().size());
        SvCluster cluster = tester.getClusters().get(0);

        assertTrue(cluster.getResolvedType() == DEL_TI);
        assertEquals(var2.position(false) - var1.position(false), getSyntheticLength(cluster));
        assertEquals(var2.position(true) - var1.position(true), getSyntheticTiLength(cluster));
        assertEquals(var1.position(false) - var2.position(true), getSyntheticGapLength(cluster));

        tester.clearClustersAndSVs();

        // test a synthetic version of the above
        var1 = createBnd(tester.nextVarId(), "1", 1000, -1, "2", 100, -1);
        var2 = createBnd(tester.nextVarId(), "1", 40000, 1, "2", 200, 1);
        SvVarData var3 = createBnd(tester.nextVarId(), "1", 36000, 1, "3", 100, -1);
        SvVarData var4 = createBnd(tester.nextVarId(), "1", 60000, -1, "3", 200, 1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);
        tester.preClusteringInit();
        tester.Analyser.clusterAndAnalyse();

        assertEquals(1, tester.Analyser.getClusters().size());
        cluster = tester.getClusters().get(0);

        assertTrue(cluster.getResolvedType() == DEL_TI);
        assertEquals(var4.position(true) - var2.position(true), getSyntheticLength(cluster));
        assertEquals(var3.position(true) - var1.position(true), getSyntheticTiLength(cluster));
        assertEquals(var2.position(true) - var3.position(true), getSyntheticGapLength(cluster));
    }

}
