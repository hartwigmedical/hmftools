package com.hartwig.hmftools.pavereverse.variants;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.hartwig.hmftools.pavereverse.base.PaddedExon;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;
import com.hartwig.hmftools.pavereverse.base.CodonWithinExons;

import org.junit.Assert;
import org.junit.Test;

public class StartLostTest extends VariantTest
{
    //MAAQVAPAAS  ->
    private final AminoAcidRange fsRange = new AminoAcidRange(aas(1, "M"), aas(1, "M"));

    @Test
    public void variantSequenceTest()
    {
        StartLost startLost = new StartLost(gene, transcript, taa, fsRange);
        Assert.assertEquals(0, startLost.variantSequence().length());
    }

    @Test
    public void applyChangeTest()
    {
        PaddedExon exon = new PaddedExon(8, "", "", exon0Bases, 9, "GGATC", "TACG");
        ChangeContext context = new ChangeContext(exon, 6, 6, true, 1);
        //MAAQVAPAAS  ->
        final AminoAcidRange range = new AminoAcidRange(aas(1, "M"), aas(1, "M"));
        StartLost startLost = new StartLost(gene, transcript, taa, range);
        // M   A   A   Q   V...
        // ATG GCC GCG CAG GTC...
        // Need ATG -> {anything else} @ 1st codon. We just report the SNVs.
        Set<ChangeResult> results = startLost.applyChange(context);
        Assert.assertEquals(9, results.size());
        AminoAcidSequence empty = AminoAcidSequence.empty();
        assertTrue(results.contains(new ChangeResult(empty, "ATAGCCGCGCAGGTC", 11, "G", "A")));
        assertTrue(results.contains(new ChangeResult(empty, "ATCGCCGCGCAGGTC", 11, "G", "C")));
        assertTrue(results.contains(new ChangeResult(empty, "ATTGCCGCGCAGGTC", 11, "G", "T")));
        assertTrue(results.contains(new ChangeResult(empty, "ACGGCCGCGCAGGTC", 10, "T", "C")));
        assertTrue(results.contains(new ChangeResult(empty, "AGGGCCGCGCAGGTC", 10, "T", "G")));
        assertTrue(results.contains(new ChangeResult(empty, "AAGGCCGCGCAGGTC", 10, "T", "A")));
        assertTrue(results.contains(new ChangeResult(empty, "CTGGCCGCGCAGGTC", 9, "A", "C")));
        assertTrue(results.contains(new ChangeResult(empty, "GTGGCCGCGCAGGTC", 9, "A", "G")));
        assertTrue(results.contains(new ChangeResult(empty, "TTGGCCGCGCAGGTC", 9, "A", "T")));
    }

    @Test
    public void possibleChangesTest()
    {
        StartLost startLost = new StartLost(gene, transcript, taa, fsRange);
        Set<CodonChange> changes = startLost.possibleVariants(new CodonWithinExons(bs(10, "ATG", true)));
        assertEquals(9, changes.size());
        assertTrue(changes.contains(new CodonChange("ATG", "ATA")));
        assertTrue(changes.contains(new CodonChange("ATG", "ATC")));
        assertTrue(changes.contains(new CodonChange("ATG", "ATT")));
        assertTrue(changes.contains(new CodonChange("ATG", "AAG")));
        assertTrue(changes.contains(new CodonChange("ATG", "ACG")));
        assertTrue(changes.contains(new CodonChange("ATG", "AGG")));
        assertTrue(changes.contains(new CodonChange("ATG", "CTG")));
        assertTrue(changes.contains(new CodonChange("ATG", "GTG")));
        assertTrue(changes.contains(new CodonChange("ATG", "TTG")));
    }
}
