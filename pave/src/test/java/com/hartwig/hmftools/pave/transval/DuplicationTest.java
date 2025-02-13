package com.hartwig.hmftools.pave.transval;

import org.junit.Assert;
import org.junit.Test;

public class DuplicationTest extends VariantTest
{
    @Test
    public void calculateVariantTest()
    {
        Duplication duplication = new Duplication(gene, transcript, taa, aar);
        TransvalVariant variant = duplication.calculateVariant(fsg);
        checkSingleHotspot(variant, "C", "CGCGCAG", "chr5", 14);
    }

    @Test
    public void variantSequenceTest()
    {
        Duplication duplication = new Duplication(gene, transcript, taa, aar);
        AminoAcidSequence expected = AminoAcidSequence.parse("MAAQAQVAPAAS");
        Assert.assertEquals(expected, duplication.variantSequence());
    }

    private Duplication dup(String gene, String variant)
    {
        return transval.variationParser().parseDuplication(gene, variant);
    }
}
