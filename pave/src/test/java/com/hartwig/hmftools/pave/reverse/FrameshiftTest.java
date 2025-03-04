package com.hartwig.hmftools.pave.reverse;

import java.util.Set;

import org.junit.Assert;
import org.junit.Test;

public class FrameshiftTest extends VariantTest
{
    //MAAQVAPAAS  -> MAAQVA*
    private final AminoAcidRange fsRange = new AminoAcidRange(aas(7, "P"), aas(7, "P"));

    @Test
    public void variantSequenceTest()
    {
        Frameshift frameshift = new Frameshift(gene, transcript, taa, fsRange);
        AminoAcidSequence expected = AminoAcidSequence.parse("MAAQVA");
        Assert.assertEquals(expected, frameshift.variantSequence());
    }

    @Test
    public void isConsistentWithThisVariantTest()
    {
        Frameshift fs = new Frameshift(gene, transcript, taa, fsRange);
        Assert.assertTrue(fs.isConsistentWithThisVariant(AminoAcidSequence.parse("MAAQVARTY")));
        Assert.assertTrue(fs.isConsistentWithThisVariant(AminoAcidSequence.parse("MAAQVATTSARTY")));
        Assert.assertTrue(fs.isConsistentWithThisVariant(AminoAcidSequence.parse("MAAQVAL")));
        // Not consistent if the P matches, as it's meant to be the first changed amino acid.
        Assert.assertFalse(fs.isConsistentWithThisVariant(AminoAcidSequence.parse("MAAQVAPKL")));
        Assert.assertFalse(fs.isConsistentWithThisVariant(AminoAcidSequence.parse("MAAQVAP")));
        // Not consistent if the start doesn't match.
        Assert.assertFalse(fs.isConsistentWithThisVariant(AminoAcidSequence.parse("MAA")));
        Assert.assertFalse(fs.isConsistentWithThisVariant(AminoAcidSequence.parse("MAAQVA"))); // Candidate needs to be longer than the variant.
        Assert.assertFalse(fs.isConsistentWithThisVariant(AminoAcidSequence.parse("MACAQVA")));
        Assert.assertFalse(fs.isConsistentWithThisVariant(AminoAcidSequence.parse("MMAAQVA")));
    }

    @Test
    public void applyChangeTest()
    {
        PaddedExon exon = new PaddedExon(8,"", "", exon0Bases, 9, "GGATC", "TACG" );
        ChangeContext context = new ChangeContext(exon, 6, 6, true, 1);
        //MAAQVAPAAS  -> MA*
        final AminoAcidRange range = new AminoAcidRange(aas(3, "A"), aas(3, "A"));
        Frameshift fs = new Frameshift(gene, transcript, taa, range);
        // M   A   A   Q   V
        // ATG GCC GCG CAG GTC
        // ATG GCC CGC AGG TCT... M A R R S...
        Set<ChangeResult> results = fs.applyChange(context);
        Assert.assertEquals(1, results.size());
        ChangeResult result = results.iterator().next();
        Assert.assertEquals("MARRS", result.mAminoAcids.sequence());
        String bases = result.mBases;
        Assert.assertEquals("ATGGCCCGCAGGTCT", bases);
        Assert.assertEquals(9 + 6 - 1, result.mLocation);
        Assert.assertEquals("CG", result.mRefBases);
        Assert.assertEquals("C", result.mAltBases);
    }
}
