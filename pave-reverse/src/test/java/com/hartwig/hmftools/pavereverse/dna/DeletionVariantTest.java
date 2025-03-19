package com.hartwig.hmftools.pavereverse.dna;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
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

public class DeletionVariantTest extends ReversePaveTestBase
{
    @Test
    public void toGenomicVariant()
    {
        int[] exonStarts = { 10, 30, 50, 70, 90 };
        int codingStart = 55;
        int codingEnd = 95;
        // 12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        // 1        10        20        30        40        50        60        70        80        90        100       110
        // ____-----++++++++++----------++++++++++----------+++++*****----------**********----------******++++------_____
        // ACTGCCCCCACGTACGTACTTTTTTTTTTACGTACGTACAAAAAAAAAACCCCCATGCCATATATATATCGTGCTAGGGTATATATATACCTTAAGGGGCCCCCCAATTC

        GeneData geneData = createEnsemblGeneData("id132", "BLAH", "1",  POS_STRAND, 5, 105);
        TranscriptData transcript = createTransExons(geneData.GeneId, 123, POS_STRAND, exonStarts, 9, codingStart, codingEnd, false, "whatever");
        RefGenomeInterface genome = new FixedStringGenome("ACTGCCCCCACGTACGTACTTTTTTTTTTACGTACGTACAAAAAAAAAACCCCCATGCCATATATATATCGGGGTAGGGTATATATATACCTTAAGGGGCCCCCCAATTC");

        DeletionVariant deletionVariant = new DeletionVariant(geneData, transcript, new InExon(7), "G");
        BaseSequenceChange change = deletionVariant.toGenomicVariant(genome);
        assertEquals("1", change.Chromosome);
        assertEquals(70, change.Position);
        assertEquals("C", change.Alt);
        assertEquals("CG", change.Ref);
    }
}
