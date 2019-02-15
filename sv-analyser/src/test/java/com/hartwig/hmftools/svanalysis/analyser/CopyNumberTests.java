package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.calcAdjustedPloidyValues;

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


}
