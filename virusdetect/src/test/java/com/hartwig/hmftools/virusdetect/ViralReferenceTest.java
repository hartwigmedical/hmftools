package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.virusdetect.ViralReference.InfoRow;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class ViralReferenceTest
{
    // Three contigs: contigA alone in its group, contigB and contigC sharing a group.
    private static final InfoRow INFO_A = new InfoRow("Virus Alpha", "Group Alpha");
    private static final InfoRow INFO_B = new InfoRow("Virus Beta type 1", "Group Beta");
    private static final InfoRow INFO_C = new InfoRow("Virus Beta type 2", "Group Beta");
    private static final Map<String, InfoRow> INFO = Map.of("contigA", INFO_A, "contigB", INFO_B, "contigC", INFO_C);

    private static final ViralContig CONTIG_A = new ViralContig("contigA", 100, "Virus Alpha", "Group Alpha");
    private static final ViralContig CONTIG_B = new ViralContig("contigB", 200, "Virus Beta type 1", "Group Beta");
    private static final ViralContig CONTIG_C = new ViralContig("contigC", 300, "Virus Beta type 2", "Group Beta");

    private static final String VALID_INFO_TSV = """
            ref_contig\tvirus_name\toncology_group
            contigA\tVirus Alpha\tGroup Alpha
            contigB\tVirus Beta type 1\tGroup Beta
            contigC\tVirus Beta type 2\tGroup Beta
            """;

    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    @Test
    public void testJoinsContigsInFastaOrderWithInfo()
    {
        List<ViralContig> contigs = ViralReference.join(dictionary("contigA", "contigB", "contigC"), INFO);
        assertEquals(List.of(CONTIG_A, CONTIG_B, CONTIG_C), contigs);
    }

    @Test
    public void testJoinThrowsWhenContigHasNoInfoRow()
    {
        Map<String, InfoRow> infoMissingC = Map.of("contigA", INFO_A, "contigB", INFO_B);
        assertThrows(UserInputError.class, () -> ViralReference.join(dictionary("contigA", "contigB", "contigC"), infoMissingC));
    }

    @Test
    public void testJoinThrowsWhenInfoRowHasNoContig()
    {
        assertThrows(UserInputError.class, () -> ViralReference.join(dictionary("contigA", "contigB"), INFO));
    }

    @Test
    public void testLoadsInfoRowsByContig() throws IOException
    {
        assertEquals(INFO, ViralReference.loadInfo(writeTsv(VALID_INFO_TSV)));
    }

    @Test
    public void testLoadInfoThrowsOnDuplicateContig() throws IOException
    {
        String duplicate = VALID_INFO_TSV + "contigA\tVirus Alpha\tGroup Alpha\n";
        assertThrows(UserInputError.class, () -> ViralReference.loadInfo(writeTsv(duplicate)));
    }

    @Test
    public void testLoadInfoThrowsOnMissingColumn() throws IOException
    {
        String noGroup = """
                ref_contig\tvirus_name
                contigA\tVirus Alpha
                """;
        assertThrows(UserInputError.class, () -> ViralReference.loadInfo(writeTsv(noGroup)));
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
