package com.hartwig.hmftools.pavereverse.dna;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.FixedStringGenome;
import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;

import org.junit.Test;

public class DnaVariantTest extends ReversePaveTestBase
{
    int[] exonStarts = { 10, 30, 50, 70, 90 };
    int codingStart = 55;
    int codingEnd = 95;
    // 12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    // 1        10        20        30        40        50        60        70        80        90        100       110
    // ____-----++++++++++----------++++++++++----------+++++*****----------**********----------******++++------_____
    // ACTGCCCCCACGTACGTACTTTTTTTTTTACGTACGTACAAAAAAAAAACCCCCATGCCATATATATATCGTGCTAGGGTATATATATACCTTAAGGGGCCCCCCAATTC

    GeneData geneData = createEnsemblGeneData("id132", "BLAH", "1", POS_STRAND, 5, 105);
    GeneData geneDataRS = createEnsemblGeneData("id133", "BLAG", "1", NEG_STRAND, 5, 105);
    TranscriptData transcript =
            createTransExons(geneData.GeneId, 123, POS_STRAND, exonStarts, 9, codingStart, codingEnd, false, "whatever");
    TranscriptData transcriptRS =
            createTransExons(geneDataRS.GeneId, 124, NEG_STRAND, exonStarts, 9, codingStart, codingEnd, false, "whatever");
    RefGenomeInterface genome =
            new FixedStringGenome("ACTGCCCCCACGTACGTACTTTTTTTTTTACGTACGTACAAAAAAAAAACCCCCATGCCATATATATATCGTGCTAGGGTATATATATACCTTAAGGGGCCCCCCAATTC");

    @Test
    public void insertionReverseStrand()
    {
        InsertionVariant duplication = new InsertionVariant(geneDataRS, transcriptRS, new InExon(8), new InExon(9), "GTA");
        BaseSequenceChange change = duplication.toGenomicVariant(genome);
        check(change,77, "G", "GTAC" );
    }

    @Test
    public void insertion()
    {
        InsertionVariant duplication = new InsertionVariant(geneData, transcript, new InExon(9), new InExon(10), "GTA");
        BaseSequenceChange change = duplication.toGenomicVariant(genome);
        check(change,73, "G", "GGTA" );

        duplication = new InsertionVariant(geneData, transcript, new InIntronAfterExon(15, 2), new InIntronAfterExon(15, 3), "GTA");
        change = duplication.toGenomicVariant(genome);
        check(change,81, "A", "AGTA" );
    }

    @Test
    public void duplicationOfRangeReverseStrand()
    {
        DuplicationVariant duplication = new DuplicationVariant(geneDataRS, transcriptRS, new InExon(9), new InExon(11), "CTA");
        BaseSequenceChange change = duplication.toGenomicVariant(genome);
        check(change,74, "C", "CTAG" );
    }

    @Test
    public void duplicationReverseStrand()
    {
        DuplicationVariant duplication = new DuplicationVariant(geneDataRS, transcriptRS, new InExon(2), new InExon(2), "T");
        BaseSequenceChange change = duplication.toGenomicVariant(genome);
        check(change,93, "T", "TA" );
    }

    @Test
    public void duplication()
    {
        DuplicationVariant duplication = new DuplicationVariant(geneData, transcript, new InExon(16), new InExon(16), "C");
        BaseSequenceChange change = duplication.toGenomicVariant(genome);
        check(change,89, "A", "AC" );
    }

    @Test
    public void duplicationOfRange()
    {
        DuplicationVariant duplication = new DuplicationVariant(geneData, transcript, new InExon(6), new InExon(10), "C");
        BaseSequenceChange change = duplication.toGenomicVariant(genome);
        check(change,69, "T", "TCGTGC" );
    }

    @Test
    public void deletionOfRange()
    {
        DeletionVariant deletionVariant = new DeletionVariant(geneData, transcript, new InExon(7), new InExon(10), "GTGC");
        BaseSequenceChange change = deletionVariant.toGenomicVariant(genome);
        check(change, 70, "CGTGC", "C");
    }

    @Test
    public void deletionOfRangeReverseStrand()
    {
        DeletionVariant deletionVariant = new DeletionVariant(geneDataRS, transcriptRS, new InExon(7), new InExon(10), "AGGG");
        BaseSequenceChange change = deletionVariant.toGenomicVariant(genome);
        check(change, 75, "TAGGG", "T");
    }

    @Test
    public void substitution()
    {
        SubstitutionVariant variant = new SubstitutionVariant(geneData, transcript, new InExon(7), "G", "T");
        BaseSequenceChange change = variant.toGenomicVariant(genome);
        check(change, 71, "G", "T");
    }

    @Test
    public void substitutionRS()
    {
        SubstitutionVariant variant = new SubstitutionVariant(geneDataRS, transcriptRS, new InExon(2), "A", "G");
        BaseSequenceChange change = variant.toGenomicVariant(genome);
        check(change, 94, "T", "C");
    }

    @Test
    public void deletion()
    {
        DeletionVariant deletionVariant = new DeletionVariant(geneData, transcript, new InExon(7), "G");
        BaseSequenceChange change = deletionVariant.toGenomicVariant(genome);
        check(change, 70, "CG", "C");
    }

    @Test
    public void deletionRS()
    {
        DeletionVariant deletionVariant = new DeletionVariant(geneDataRS, transcriptRS, new InExon(14), "A");
        BaseSequenceChange change = deletionVariant.toGenomicVariant(genome);
        check(change, 71, "GT", "G");
    }

    private static void check(final BaseSequenceChange change, final int expected, final String ref, final String alt)
    {
        assertEquals("1", change.Chromosome);
        assertEquals(expected, change.Position);
        assertEquals(ref, change.Ref);
        assertEquals(alt, change.Alt);
    }
}
