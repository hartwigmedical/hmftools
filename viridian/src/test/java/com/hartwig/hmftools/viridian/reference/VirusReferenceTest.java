package com.hartwig.hmftools.viridian.reference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.viridian.common.UserInputError;

import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class VirusReferenceTest
{
    private static final VirusGroup GROUP_ALPHA = new VirusGroup("Group Alpha");
    private static final VirusGroup GROUP_BETA = new VirusGroup("Group Beta");

    // Three virus genomes: contigA alone in its group, contigB and contigC sharing a group.
    private static final VirusGenomeInfo INFO_A =
            new VirusGenomeInfo("contigA", new Virus("Virus Alpha", GROUP_ALPHA));
    private static final VirusGenomeInfo INFO_B =
            new VirusGenomeInfo("contigB", new Virus("Virus Beta type 1", GROUP_BETA));
    private static final VirusGenomeInfo INFO_C =
            new VirusGenomeInfo("contigC", new Virus("Virus Beta type 2", GROUP_BETA));
    private static final List<VirusGenomeInfo> INFO = List.of(INFO_A, INFO_B, INFO_C);

    private static final VirusGenome GENOME_A = new VirusGenome("contigA", 100, INFO_A.virus());
    private static final VirusGenome GENOME_B = new VirusGenome("contigB", 200, INFO_B.virus());
    private static final VirusGenome GENOME_C = new VirusGenome("contigC", 300, INFO_C.virus());

    @Test
    public void testJoinFastaAndInfoJoinsGenomesInFastaOrder()
    {
        List<VirusGenome> genomes = VirusReference.joinFastaAndInfo(dictionary(GENOME_A, GENOME_B, GENOME_C), INFO);
        assertEquals(List.of(GENOME_A, GENOME_B, GENOME_C), genomes);
    }

    @Test
    public void testJoinFastaAndInfoThrowsWhenContigHasNoInfoRow()
    {
        SAMSequenceDictionary dictionary = dictionary(GENOME_A, GENOME_B, GENOME_C);
        List<VirusGenomeInfo> infoMissingC = List.of(INFO_A, INFO_B);

        assertThrows(UserInputError.class, () -> VirusReference.joinFastaAndInfo(dictionary, infoMissingC));
    }

    @Test
    public void testJoinFastaAndInfoThrowsWhenInfoRowHasNoContig()
    {
        SAMSequenceDictionary dictionaryMissingC = dictionary(GENOME_A, GENOME_B);

        assertThrows(UserInputError.class, () -> VirusReference.joinFastaAndInfo(dictionaryMissingC, INFO));
    }

    private static SAMSequenceDictionary dictionary(VirusGenome... genomes)
    {
        return new SAMSequenceDictionary(Arrays.stream(genomes)
                .map(genome -> new SAMSequenceRecord(genome.contigName(), genome.length()))
                .toList());
    }
}
