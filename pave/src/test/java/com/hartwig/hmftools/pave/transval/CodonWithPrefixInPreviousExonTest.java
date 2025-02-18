package com.hartwig.hmftools.pave.transval;

import org.junit.Assert;
import org.junit.Test;

public class CodonWithPrefixInPreviousExonTest extends TransvalTest
{
    @Test
    public void locationTest()
    {
        CodonWithPrefixInPreviousExon codon = new CodonWithPrefixInPreviousExon("A", bs(100, "TC"));
        Assert.assertEquals(100, codon.strandLocationOfChange(1));
        Assert.assertEquals(101, codon.strandLocationOfChange(2));

        codon = new CodonWithPrefixInPreviousExon("AG", bs(100, "T"));
        Assert.assertEquals(100, codon.strandLocationOfChange(2));
    }
}
