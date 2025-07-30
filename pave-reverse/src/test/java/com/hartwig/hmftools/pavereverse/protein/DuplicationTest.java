package com.hartwig.hmftools.pavereverse.protein;

import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;

import org.junit.Assert;
import org.junit.Test;

public class DuplicationTest extends VariantTest
{
    @Test
    public void calculateVariantTest()
    {
        Duplication duplication = new Duplication(gene, transcript, taa, aar);
        BaseSequenceVariants variant = duplication.calculateVariant(fsg);
        checkSingleChange(variant, "C", "CGCGCAG", "chr5", 15);
    }

    @Test
    public void variantSequenceTest()
    {
        Duplication duplication = new Duplication(gene, transcript, taa, aar);
        AminoAcidSequence expected = AminoAcidSequence.parse("MAAQAQVAPAAS");
        Assert.assertEquals(expected, duplication.variantSequence());
    }
}
