package com.hartwig.hmftools.sig_analyser;

import static com.hartwig.hmftools.sig_analyser.common.DataUtils.scaleRoundRatio;

public class MiscTest
{

    private static void testRoundFunction()
    {
        double value = 0.25214;
        double result = scaleRoundRatio(value, 1);
        result = scaleRoundRatio(value, 2);
        result = scaleRoundRatio(value, 3);

        value = 0.75214; // should scale to be rounded to 0.1
        result = scaleRoundRatio(value, 1);
        result = scaleRoundRatio(value, 2);

    }

}
