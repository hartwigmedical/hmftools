package com.hartwig.hmftools.pave.transval;

import org.junit.Assert;
import org.junit.Test;

public class CodonWithSuffixInNextExonTest extends TransvalTest
{
    @Test
    public void locationTest()
    {
        CodonWithSuffixInNextExon codon = new CodonWithSuffixInNextExon(bs(100, "TC", true), "G");
        Assert.assertEquals(100, codon.forwardStrandLocationOfChange(0));
        Assert.assertEquals(101, codon.forwardStrandLocationOfChange(1));
    }
}
