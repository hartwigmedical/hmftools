package com.hartwig.hmftools.pavereverse.variants;

import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.hartwig.hmftools.pavereverse.base.PaddedExon;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;

import org.junit.Assert;
import org.junit.Test;

public class StopGainedTest extends VariantTest
{
    //MAAQVAPAAS  -> MAA
    private final AminoAcidRange fsRange = new AminoAcidRange(aas(4, "Q"), aas(4, "Q"));

    @Test
    public void variantSequenceTest()
    {
        StopGained stopGained = new StopGained(gene, transcript, taa, fsRange);
        AminoAcidSequence expected = AminoAcidSequence.parse("MAA");
        Assert.assertEquals(expected, stopGained.variantSequence());
    }

    @Test
    public void applyChangeTest()
    {
        PaddedExon exon = new PaddedExon(8, "", "", exon0Bases, 9, "GGATC", "TACG");
        ChangeContext context = new ChangeContext(exon, 6, 6, true, 1);
        //MAAQVAPAAS  -> MAA
        final AminoAcidRange range = new AminoAcidRange(aas(4, "Q"), aas(4, "Q"));
        StopGained stopGained = new StopGained(gene, transcript, taa, range);
        // M   A   A   Q   V...
        // ATG GCC GCG CAG GTC...
        // Need CAG -> {TAA, TAG, TGA} @ fourth codon
        Set<ChangeResult> results = stopGained.applyChange(context);
        Assert.assertEquals(3, results.size());
        assertTrue(results.contains(new ChangeResult(aaSeq("MAAX"), "ATGGCCGCGTAGGTC", 18, "C", "T")));
        assertTrue(results.contains(new ChangeResult(aaSeq("MAAX"), "ATGGCCGCGTGAGTC", 18, "CAG", "TGA")));
        assertTrue(results.contains(new ChangeResult(aaSeq("MAAX"), "ATGGCCGCGTAAGTC", 18, "CAG", "TAA")));
    }
}
