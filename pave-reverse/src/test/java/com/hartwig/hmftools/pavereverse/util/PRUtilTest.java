package com.hartwig.hmftools.pavereverse.util;

import org.junit.Assert;
import org.junit.Test;

public class PRUtilTest
{
    @Test
    public void substitutionDistanceTest()
    {
        Assert.assertEquals(0, PRUtils.substitutionDistance("", ""));
        Assert.assertEquals(0, PRUtils.substitutionDistance("A", "A"));
        Assert.assertEquals(0, PRUtils.substitutionDistance("AA", "AA"));
        Assert.assertEquals(1, PRUtils.substitutionDistance("A", "B"));
        Assert.assertEquals(3, PRUtils.substitutionDistance("ABCDEAB", "ABGDTTB"));
    }
}
