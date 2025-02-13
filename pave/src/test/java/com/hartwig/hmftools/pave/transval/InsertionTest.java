package com.hartwig.hmftools.pave.transval;

import java.util.Set;

import com.hartwig.hmftools.common.codon.Nucleotides;

import org.junit.Assert;
import org.junit.Test;

public class InsertionTest extends VariantTest
{
    @Test
    public void calculateVariantTest()
    {
//        Duplication duplication = new Duplication(gene, transcript, taa, aar);
//        TransvalVariant variant = duplication.calculateVariant(fsg);
//        checkSingleHotspot(variant, "C", "CGCGCAG", "chr5", 14);
    }

//    @Test
    public void applyChangeTest()
    {
        PaddedExon exon = new PaddedExon("", "", exon0Bases, 9, "GGATC" );
        ChangeContext context = new ChangeContext(exon, 3, 4, true, 1);
        Insertion insertion = new Insertion(gene, transcript, taa, aar, aaSeq("SPQR"));
        ChangeResult result = insertion.applyChange(context);
        Assert.assertEquals("MAASPQRQVAPAAS", result.mAminoAcids.sequence());
        String bases = result.mBases;
        Assert.assertTrue(bases.startsWith(exon0Bases.substring(0,3)));
        Assert.assertTrue(bases.endsWith(exon0Bases.substring(4)));
        Assert.assertEquals("SPQR", AminoAcidSequence.fromNucleotides(bases.substring(3, 14)).sequence());
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
