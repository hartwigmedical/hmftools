package com.hartwig.hmftools.purple.fitting;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedPurity;

import org.junit.Test;

public class BestFitTest extends FittingTestBase
{
    @Test
    public void testMostDiploidPurity()
    {
        FittedPurity fp1 = createFittedPurity(0.3, 0.3, 2.3);
        FittedPurity fp2 = createFittedPurity(0.3, 0.2, 1.9);
        FittedPurity fp3 = createFittedPurity(0.4, 0.4, 1.8);
        FittedPurity fp4 = createFittedPurity(0.4, 0.3, 2.05);

        List<FittedPurity> all = Lists.newArrayList(fp1, fp2, fp3, fp4);
        Collections.shuffle(all);

        List<FittedPurity> result = BestFit.mostDiploidPerPurity(all);

        assertEquals(2, result.size());
        assertEquals(fp2, result.get(0));
        assertEquals(fp4, result.get(1));
    }
}