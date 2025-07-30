package com.hartwig.hmftools.pavereverse.base;

import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;

import org.junit.Assert;
import org.junit.Test;

public class CodonWithSuffixInNextExonTest extends ReversePaveTestBase
{
    @Test
    public void locationTest()
    {
        CodonWithSuffixInNextExon codon = new CodonWithSuffixInNextExon(bs(100, "TC", true), "G");
        Assert.assertEquals(100, codon.forwardStrandLocationOfChange(0));
        Assert.assertEquals(101, codon.forwardStrandLocationOfChange(1));
    }
}
