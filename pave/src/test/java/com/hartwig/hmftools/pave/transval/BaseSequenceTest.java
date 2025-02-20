package com.hartwig.hmftools.pave.transval;

import org.junit.Assert;
import org.junit.Test;

public class BaseSequenceTest extends TransvalTestBase
{
    @Test
    public void reverseComplementTest()
    {
        BaseSequence bs = new BaseSequence(100, "A", true);
        Assert.assertEquals(100, bs.reverseComplement().Start);
        Assert.assertEquals("T", bs.reverseComplement().Bases);
        Assert.assertFalse(bs.reverseComplement().IsForwardStrand);
        Assert.assertEquals(bs, bs.reverseComplement().reverseComplement());

        bs = new BaseSequence(100, "AC", true);
        Assert.assertEquals(101, bs.reverseComplement().Start);
        Assert.assertEquals("GT", bs.reverseComplement().Bases);
        Assert.assertFalse(bs.reverseComplement().IsForwardStrand);
        Assert.assertEquals(bs, bs.reverseComplement().reverseComplement());

        bs = new BaseSequence(100, "ACG", true);
        Assert.assertEquals(102, bs.reverseComplement().Start);
        Assert.assertEquals("CGT", bs.reverseComplement().Bases);
        Assert.assertFalse(bs.reverseComplement().IsForwardStrand);
        Assert.assertEquals(bs, bs.reverseComplement().reverseComplement());
    }
}
