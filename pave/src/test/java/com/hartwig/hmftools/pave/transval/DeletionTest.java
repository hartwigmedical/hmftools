package com.hartwig.hmftools.pave.transval;

import org.junit.Assert;
import org.junit.Test;

public class DeletionTest extends VariantTest
{
    @Test
    public void variantSequenceTest()
    {
        Deletion deletion = new Deletion(gene, transcript, taa, aar);
        AminoAcidSequence expected = AminoAcidSequence.parse("MAVAPAAS");
        Assert.assertEquals(expected, deletion.variantSequence());
    }
}
