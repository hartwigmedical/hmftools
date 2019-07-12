package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.linx.analyser.SvTestHelper.createDel;
import static com.hartwig.hmftools.linx.cn.CnPloidyCalcs.calcAdjustedPloidyValues;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import com.hartwig.hmftools.linx.cn.PloidyCalcData;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

public class CopyNumberTest
{
    @Test
    public void testPloidyRecalcs()
    {
        // example
        double cnChgStart = 3.1;
        double cnChgEnd = 3.6;
        int tumorReadCount = 58;
        double maxCNStart = 6.01;
        double maxCNEnd = 6.66;
        double ploidy = 2.14;
        int[] startDepthData = {3, 3};
        int[] endDepthData = {1, 1};

        PloidyCalcData calcResults =  calcAdjustedPloidyValues(cnChgStart, cnChgEnd, tumorReadCount, ploidy,
                maxCNStart, maxCNEnd, startDepthData, endDepthData);

        assertEquals(2.88, calcResults.PloidyEstimate, 0.01);
        assertEquals(0.80, calcResults.PloidyUncertainty, 0.01);

        cnChgStart = 2;
        cnChgEnd = 1.3;
        tumorReadCount = 69;
        maxCNStart = 7;
        maxCNEnd = 5;
        ploidy = 1.43;
        startDepthData[0] = 0;;
        endDepthData[0] = 0;

        calcResults =  calcAdjustedPloidyValues(cnChgStart, cnChgEnd, tumorReadCount, ploidy,
                maxCNStart, maxCNEnd, startDepthData, endDepthData);

        assertEquals(1.49, calcResults.PloidyEstimate, 0.01);
        assertEquals(0.79, calcResults.PloidyUncertainty, 0.01);
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

        assertFalse(cluster1.hasVariedPloidy());
        assertEquals(1.0, cluster1.getMinPloidy(), 0.001);
        assertEquals(1.0, cluster1.getMaxPloidy(), 0.001);

        // now check alignment to a common but non-integer ploidy
        SvVarData var4 = createDel(tester.nextVarId(), "1", 700, 800);

        // all SVs have 1.1-1.2 in common
        var4.setPloidyRecalcData(1.1, 1.5);
        cluster1.addVariant(var4);

        assertFalse(cluster1.hasVariedPloidy());
        assertEquals(1.0, cluster1.getMinPloidy(), 0.001);
        assertEquals(1.0, cluster1.getMaxPloidy(), 0.001);


    }



}
