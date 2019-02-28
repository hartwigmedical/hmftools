package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.calcAdjustedPloidyValues;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Test;

public class CopyNumberTests
{
    @Test
    public void testPloidyRecalcs()
    {
        // example
        double cnStart = 3;
        double cnChgStart = 0;
        double cnEnd = 3;
        double cnChgEnd = 4;
        double ploidy = 4;
        int tumorReadCount = 25;
        int[] startDepthData = {100, 100};
        int[] endDepthData = {1000, 1000};

        double[] calcResults =  calcAdjustedPloidyValues(cnStart, cnChgStart, cnEnd, cnChgEnd, ploidy,
                tumorReadCount, startDepthData, endDepthData);

        // assertEqual()

        // prod example
        cnStart = 5.684;
        cnChgStart = 0.941;
        cnEnd = 5.764;
        cnChgEnd = 1.021;
        ploidy = 1.017;
        tumorReadCount = 22;
        endDepthData[0] = 3;
        endDepthData[1] = 0;

        startDepthData[0] = 0;
        startDepthData[1] = 0;

        calcResults =  calcAdjustedPloidyValues(cnStart, cnChgStart, cnEnd, cnChgEnd, ploidy,
                tumorReadCount, startDepthData, endDepthData);


    }

    @Test
    public void testClusterMinMaxCopyNumber()
    {
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        SvVarData var1 = createDel(tester.nextVarId(), "1", 100, 200);
        SvVarData var2 = createDel(tester.nextVarId(), "1", 300, 400);
        SvVarData var3 = createDel(tester.nextVarId(), "1", 500, 600);

        var1.setPloidyRecalcData(0.9, 1.6);
        var2.setPloidyRecalcData(0.1, 1.2);
        var3.setPloidyRecalcData(0.99, 2.01);

        SvCluster cluster1 = new SvCluster(0);
        cluster1.addVariant(var1);
        cluster1.addVariant(var2);
        cluster1.addVariant(var3);

        assertFalse(cluster1.hasVariedCopyNumber());
        assertEquals(1.0, cluster1.getMinCNChange(), 0.001);
        assertEquals(1.0, cluster1.getMaxCNChange(), 0.001);

        // now check alignment to a common but non-integer ploidy
        SvVarData var4 = createDel(tester.nextVarId(), "1", 700, 800);

        // all SVs have 1.1-1.2 in common
        var4.setPloidyRecalcData(1.1, 1.5);
        cluster1.addVariant(var4);

        assertFalse(cluster1.hasVariedCopyNumber());
        assertEquals(1.0, cluster1.getMinCNChange(), 0.001);
        assertEquals(1.0, cluster1.getMaxCNChange(), 0.001);


    }



}
