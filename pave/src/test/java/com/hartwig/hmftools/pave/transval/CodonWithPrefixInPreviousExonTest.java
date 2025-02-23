package com.hartwig.hmftools.pave.transval;

import org.junit.Assert;
import org.junit.Test;

public class CodonWithPrefixInPreviousExonTest extends TransvalTestBase
{
    @Test
    public void locationTest()
    {
        CodonWithPrefixInPreviousExon codon = new CodonWithPrefixInPreviousExon("A", bs(100, "TC", true));
        Assert.assertEquals(100, codon.forwardStrandLocationOfChange(1));
        Assert.assertEquals(101, codon.forwardStrandLocationOfChange(2));

        codon = new CodonWithPrefixInPreviousExon("AG", bs(100, "T", true));
        Assert.assertEquals(100, codon.forwardStrandLocationOfChange(2));
    }
}
