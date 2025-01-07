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
        SortedSet<CodonVariant> variants = calculator.possibleCodonsForVariant("GAC");
        Assert.assertEquals(2, variants.size());
        Iterator<CodonVariant> iterator = variants.iterator();
        Assert.assertEquals(new CodonVariant("GAC", "GAA"), iterator.next());
        Assert.assertEquals(new CodonVariant("GAC", "GAG"), iterator.next());
    }
}
