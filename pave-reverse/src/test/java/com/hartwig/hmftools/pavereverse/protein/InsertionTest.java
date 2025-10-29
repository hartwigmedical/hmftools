package com.hartwig.hmftools.pavereverse.protein;

import java.util.Iterator;
import java.util.Set;

import com.hartwig.hmftools.pavereverse.base.PaddedExon;
import com.hartwig.hmftools.pavereverse.aminoacids.AminoAcidSequence;
import com.hartwig.hmftools.pavereverse.base.ChangeContext;

import org.junit.Assert;
import org.junit.Test;

public class InsertionTest extends VariantTest
{
    @Test
    public void applySimpleChangeTest()
    {
        PaddedExon exon = new PaddedExon(8,"", "", exon0Bases, 9, "GGATC", "TACG" );
        ChangeContext context = new ChangeContext(exon, 9, 10, true, 1);
        Insertion insertion = new Insertion(gene, transcript, taa, aar, aaSeq("Y"));
        Set<ChangeResult> results = insertion.applyChange(context);
        Assert.assertEquals(2, results.size());

        final Iterator<ChangeResult> iterator = results.iterator();
        ChangeResult result1 = iterator.next();
        Assert.assertEquals("MAAYQV", result1.Acids.sequence());
        String bases1 = result1.ExonBases;
        Assert.assertTrue(bases1.startsWith(exon0Bases.substring(0,9)));
        Assert.assertTrue(bases1.endsWith(exon0Bases.substring(9)));
        Assert.assertEquals("Y", AminoAcidSequence.fromNucleotides(bases1.substring(9, 12)).sequence());

        ChangeResult result2 = iterator.next();
        Assert.assertEquals("MAAYQV", result2.Acids.sequence());
        String bases2 = result2.ExonBases;
        Assert.assertTrue(bases2.startsWith(exon0Bases.substring(0,9)));
        Assert.assertTrue(bases2.endsWith(exon0Bases.substring(9)));
        Assert.assertEquals("Y", AminoAcidSequence.fromNucleotides(bases2.substring(9, 12)).sequence());

        Assert.assertNotEquals(bases1, bases2);
    }

    @Test
    public void applyChangeTest()
    {
        PaddedExon exon = new PaddedExon(8,"", "", exon0Bases, 9, "GGATC", "TACG");
        ChangeContext context = new ChangeContext(exon, 9, 10, true, 1);
        Insertion insertion = new Insertion(gene, transcript, taa, aar, aaSeq("SPQR"));
        Set<ChangeResult> results = insertion.applyChange(context);
        Assert.assertEquals(1, results.size());
        ChangeResult result = results.iterator().next();
        Assert.assertEquals("MAASPQRQV", result.Acids.sequence());
        String bases = result.ExonBases;
        Assert.assertTrue(bases.startsWith(exon0Bases.substring(0,9)));
        Assert.assertTrue(bases.endsWith(exon0Bases.substring(9)));
        Assert.assertEquals("SPQR", AminoAcidSequence.fromNucleotides(bases.substring(9, 21)).sequence());
    }

    @Test
    public void possibleInsertedNucleotideSequencesTest()
    {
        Insertion insertion = new Insertion(gene, transcript, taa, aar, aaSeq("F"));
        Set<String> expected = Set.of("TTT", "TTC");
        Assert.assertEquals(expected, insertion.possibleInsertedNucleotideSequences());

        insertion = new Insertion(gene, transcript, taa, aar, aaSeq("R"));
        expected = Set.of("AGG", "AGA", "CGA", "CGC", "CGG", "CGT");
        Assert.assertEquals(expected, insertion.possibleInsertedNucleotideSequences());

        // With more than 1 AA inserted, we just choose any single sequence that produces it.
        insertion = new Insertion(gene, transcript, taa, aar, aaSeq("FC"));
        expected = Set.of("TTTTGC", "TTTTGT", "TTCTGC", "TTCTGT");
        Set<String> actual = insertion.possibleInsertedNucleotideSequences();
        Assert.assertEquals(1, actual.size());
        Assert.assertTrue(expected.contains(actual.iterator().next()));
    }

    @Test
    public void variantSequenceTest()
    {
        Insertion insertion = new Insertion(gene, transcript, taa, aar, aaSeq("GLAM"));
        Assert.assertEquals("MAAGLAMQVAPAAS", insertion.variantSequence().sequence());
    }

    @Test
    public void positionOfFirstAlteredCodonTest()
    {
        Insertion insertion = new Insertion(gene, transcript, taa, aar, aaSeq("GLAM"));
        Assert.assertEquals(4, insertion.positionOfFirstAlteredCodon());
    }
}
