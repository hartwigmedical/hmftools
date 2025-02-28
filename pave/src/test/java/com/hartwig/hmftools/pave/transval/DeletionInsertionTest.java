package com.hartwig.hmftools.pave.transval;

import org.junit.Assert;
import org.junit.Test;

public class DeletionInsertionTest extends VariantTest
{
    @Test
    public void variantSequenceTest()
    {
        AminoAcidSequence replacement = AminoAcidSequence.parse("TKWRF");
        DeletionInsertion di = new DeletionInsertion(gene, transcript, taa, aar, replacement);
        // MAAQVAPAAS -> MA + TKWRF + VAPAAS
        AminoAcidSequence expected = AminoAcidSequence.parse("MATKWRFVAPAAS");
        Assert.assertEquals(expected, di.variantSequence());
    }

    @Test
    public void seekResultsInCompanionContextTest()
    {
        AminoAcidSequence replacement = AminoAcidSequence.parse("TKWRF");
        DeletionInsertion di = new DeletionInsertion(gene, transcript, taa, aar, replacement);
        Assert.assertFalse(di.seekResultsInCompanionContext(true));
    }
}
