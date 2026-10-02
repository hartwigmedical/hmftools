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
    private static final OncologyGroup GROUP_ALPHA = new OncologyGroup("Group Alpha");
    private static final OncologyGroup GROUP_BETA = new OncologyGroup("Group Beta");

    // Three contigs: contigA alone in its group, contigB and contigC sharing a group.
    private static final VirusInfo INFO_A = new VirusInfo("contigA", "Virus Alpha", GROUP_ALPHA);
    private static final VirusInfo INFO_B = new VirusInfo("contigB", "Virus Beta type 1", GROUP_BETA);
    private static final VirusInfo INFO_C = new VirusInfo("contigC", "Virus Beta type 2", GROUP_BETA);
    private static final List<VirusInfo> INFO = List.of(INFO_A, INFO_B, INFO_C);

    private static final ViralContig CONTIG_A = new ViralContig("contigA", 100, "Virus Alpha", GROUP_ALPHA);
    private static final ViralContig CONTIG_B = new ViralContig("contigB", 200, "Virus Beta type 1", GROUP_BETA);
    private static final ViralContig CONTIG_C = new ViralContig("contigC", 300, "Virus Beta type 2", GROUP_BETA);

    private static final String VALID_INFO_TSV = """
            ref_contig\tvirus_name\toncology_group
            contigA\tVirus Alpha\tGroup Alpha
            contigB\tVirus Beta type 1\tGroup Beta
            contigC\tVirus Beta type 2\tGroup Beta
            """;

    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    @Test
    public void testJoinFastaAndInfoJoinsContigsInFastaOrder()
    {
        List<ViralContig> contigs = VirusReference.joinFastaAndInfo(dictionary("contigA", "contigB", "contigC"), INFO);
        assertEquals(List.of(CONTIG_A, CONTIG_B, CONTIG_C), contigs);
    }

    @Test
    public void testJoinFastaAndInfoThrowsWhenContigHasNoInfoRow()
    {
        List<VirusInfo> infoMissingC = List.of(INFO_A, INFO_B);
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
        assertEquals(INFO, VirusInfo.load(writeTsv(VALID_INFO_TSV)));
    }

    @Test
    public void testLoadInfoThrowsOnDuplicateContig()
    {
        String duplicate = VALID_INFO_TSV + "contigA\tVirus Alpha\tGroup Alpha\n";
        assertThrows(UserInputError.class, () -> VirusInfo.load(writeTsv(duplicate)));
    }

    @Test
    public void testLoadInfoThrowsOnMissingColumn()
    {
        String noGroup = """
                ref_contig\tvirus_name
                contigA\tVirus Alpha
                """;
        assertThrows(RuntimeException.class, () -> VirusInfo.load(writeTsv(noGroup)));
    }

    // Contig lengths are fixed per name so joined ViralContig values are predictable.
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
