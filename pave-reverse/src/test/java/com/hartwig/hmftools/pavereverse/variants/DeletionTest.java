package com.hartwig.hmftools.pavereverse.variants;

import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;

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

    @Test
    public void isConsistentWithThisVariantTest()
    {
        Deletion deletion = new Deletion(gene, transcript, taa, aar);
        Assert.assertTrue(deletion.isConsistentWithThisVariant(AminoAcidSequence.parse("MAVAPAAS")));
        Assert.assertFalse(deletion.isConsistentWithThisVariant(AminoAcidSequence.parse("MAAPAAS")));
        Assert.assertFalse(deletion.isConsistentWithThisVariant(AminoAcidSequence.parse("MAVAPAA")));
        Assert.assertFalse(deletion.isConsistentWithThisVariant(AminoAcidSequence.parse("AVAPAAS")));
    }
}
