package com.hartwig.hmftools.viridian.reference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.viridian.common.UserInputError;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class VirusGenomeInfoTest
{
    private static final String VALID_INFO_TSV = """
            contig_name\tvirus_name\tvirus_group
            contigA\tVirus Alpha\tGroup Alpha
            contigB\tVirus Beta type 1\tGroup Beta
            contigC\tVirus Beta type 2\tGroup Beta
            """;

    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    @Test
    public void testLoad() throws IOException
    {
        VirusGroup groupAlpha = new VirusGroup("Group Alpha");
        VirusGroup groupBeta = new VirusGroup("Group Beta");
        List<VirusGenomeInfo> expected = List.of(
                new VirusGenomeInfo("contigA", new Virus("Virus Alpha", groupAlpha)),
                new VirusGenomeInfo("contigB", new Virus("Virus Beta type 1", groupBeta)),
                new VirusGenomeInfo("contigC", new Virus("Virus Beta type 2", groupBeta)));

        assertEquals(expected, VirusGenomeInfo.load(writeTsv(VALID_INFO_TSV)));
    }

    @Test
    public void testLoadDuplicateContig()
    {
        String duplicate = VALID_INFO_TSV + "contigA\tVirus Alpha\tGroup Alpha\n";
        assertThrows(UserInputError.class, () -> VirusGenomeInfo.load(writeTsv(duplicate)));
    }

    @Test
    public void testLoadVirusInTwoGroups()
    {
        String virusInTwoGroups = VALID_INFO_TSV + "contigD\tVirus Alpha\tGroup Beta\n";
        assertThrows(UserInputError.class, () -> VirusGenomeInfo.load(writeTsv(virusInTwoGroups)));
    }

    private String writeTsv(String content) throws IOException
    {
        File file = new File(mTempDir.newFolder(), "info.tsv");
        Files.writeString(file.toPath(), content);
        return file.getPath();
    }
}
