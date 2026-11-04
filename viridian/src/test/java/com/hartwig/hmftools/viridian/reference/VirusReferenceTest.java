package com.hartwig.hmftools.viridian.reference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.common.UserInputError;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

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

    private static final VirusGenome GENOME_A = new VirusGenome("contigA", 100, new Virus("Virus Alpha", GROUP_ALPHA));
    private static final VirusGenome GENOME_B =
            new VirusGenome("contigB", 200, new Virus("Virus Beta type 1", GROUP_BETA));
    private static final VirusGenome GENOME_C =
            new VirusGenome("contigC", 300, new Virus("Virus Beta type 2", GROUP_BETA));

    private static final String VALID_INFO_TSV = """
            contig_name\tvirus_name\tvirus_group
            contigA\tVirus Alpha\tGroup Alpha
            contigB\tVirus Beta type 1\tGroup Beta
            contigC\tVirus Beta type 2\tGroup Beta
            """;

    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    @Test
    public void testJoinFastaAndInfoJoinsGenomesInFastaOrder()
    {
        List<VirusGenome> genomes = VirusReference.joinFastaAndInfo(dictionary("contigA", "contigB", "contigC"), INFO);
        assertEquals(List.of(GENOME_A, GENOME_B, GENOME_C), genomes);
    }

    @Test
    public void testJoinFastaAndInfoThrowsWhenContigHasNoInfoRow()
    {
        List<VirusGenomeInfo> infoMissingC = List.of(INFO_A, INFO_B);
        assertThrows(
                UserInputError.class, () -> VirusReference.joinFastaAndInfo(dictionary("contigA", "contigB", "contigC"), infoMissingC));
    }

    @Test
    public void testJoinFastaAndInfoThrowsWhenInfoRowHasNoContig()
    {
        assertThrows(UserInputError.class, () -> VirusReference.joinFastaAndInfo(dictionary("contigA", "contigB"), INFO));
    }

    @Test
    public void testLoadInfoRowsByContig() throws IOException
    {
        assertEquals(INFO, VirusGenomeInfo.load(writeTsv(VALID_INFO_TSV)));
    }

    @Test
    public void testLoadInfoThrowsOnDuplicateContig()
    {
        String duplicate = VALID_INFO_TSV + "contigA\tVirus Alpha\tGroup Alpha\n";
        assertThrows(UserInputError.class, () -> VirusGenomeInfo.load(writeTsv(duplicate)));
    }

    @Test
    public void testLoadInfoThrowsWhenVirusIsInTwoGroups()
    {
        String virusInTwoGroups = VALID_INFO_TSV + "contigD\tVirus Alpha\tGroup Beta\n";
        assertThrows(UserInputError.class, () -> VirusGenomeInfo.load(writeTsv(virusInTwoGroups)));
    }

    @Test
    public void testLoadInfoThrowsOnMissingColumn()
    {
        String noGroup = """
                contig_name\tvirus_name
                contigA\tVirus Alpha
                """;
        assertThrows(RuntimeException.class, () -> VirusGenomeInfo.load(writeTsv(noGroup)));
    }

    // Contig lengths are fixed per name so joined VirusGenome values are predictable.
    private static SAMSequenceDictionary dictionary(String... contigs)
    {
        Map<String, Integer> lengths = Map.of("contigA", 100, "contigB", 200, "contigC", 300);
        SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        for(String contig : contigs)
        {
            dictionary.addSequence(new SAMSequenceRecord(contig, lengths.get(contig)));
        }
        return dictionary;
    }

    private String writeTsv(String content) throws IOException
    {
        File file = new File(mTempDir.newFolder(), "info.tsv");
        Files.writeString(file.toPath(), content);
        return file.getPath();
    }
}
