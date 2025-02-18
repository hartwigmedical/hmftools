package com.hartwig.hmftools.pave.transval;

import java.util.Iterator;
import java.util.SortedSet;

import org.junit.Assert;
import org.junit.Test;

public class CodonVariantsCalculatorTest
{
    @Test
    public void variantsTest()
    {
        CodonVariantsCalculator calculator = new CodonVariantsCalculator("D", "E");
        SortedSet<CodonChange> variants = calculator.possibleCodonsForVariant("GAC");
        Assert.assertEquals(2, variants.size());
        Iterator<CodonChange> iterator = variants.iterator();
        Assert.assertEquals(new CodonChange("GAC", "GAA"), iterator.next());
        Assert.assertEquals(new CodonChange("GAC", "GAG"), iterator.next());
    }
}
