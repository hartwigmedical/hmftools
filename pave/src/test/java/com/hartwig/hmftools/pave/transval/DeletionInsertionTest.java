package com.hartwig.hmftools.pave.transval;

import org.junit.Assert;
import org.junit.Test;

public class DeletionInsertionTest extends TransvalTestBase
{

    @Test
    public void referenceAminoAcidsTest()
    {
        Assert.assertEquals("M", this.di("ADCK2:p.M1_M1delinsKQ").referenceAminoAcids());
        Assert.assertEquals("MVAP", this.di("ADCK2:p.M1_P4delinsK").referenceAminoAcids());
        Assert.assertEquals("EAT", this.di("ADCK2:p.E301_T303delinsQQ").referenceAminoAcids());
        Assert.assertEquals("PX", this.di("ADCK2:p.E626_T627delinsQQ").referenceAminoAcids());
    }

    @Test
    public void altAminoAcidsTest()
    {
        Assert.assertEquals("KQ", this.di("ADCK2:p.M1_M1delinsKQ").altAminoAcidSequence());
        Assert.assertEquals("K", this.di("ADCK2:p.M1_P4delinsK").altAminoAcidSequence());
    }

    private DeletionInsertion di(String definition)
    {
        return transval.variationParser().parseDeletionInsertion(definition);
    }
}
