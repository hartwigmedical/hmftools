package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.purple.tools.GermlineGeneAnalyser.MIN_FREQUENCY_CFG;
import static com.hartwig.hmftools.purple.tools.GermlineGeneAnalyser.PURPLE_GERMLINE_GENE_DATA_CSV;

import java.io.File;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class GermlineGeneAnalyserTest extends ToolsTestBase
{
    @Test
    public void mainTest() throws Exception
    {
        File samplesFile = new File(mPurpleDataDir, "samples.csv");
        File tempDir = Files.createTempDirectory("purple-test").toFile();
        tempDir.deleteOnExit();

        String[] args = new String[] {
                "-" + SAMPLE_ID_FILE, samplesFile.getAbsolutePath(),
                "-" + PURPLE_DIR_CFG, mPurpleDataDir,
                "-" + OUTPUT_DIR, tempDir.getAbsolutePath(),
                "-" + REF_GENOME_VERSION, "38",
                "-" + MIN_FREQUENCY_CFG, "2",
                "-threads", "1" };
        GermlineGeneAnalyser.main(args);

        File outputFile = new File(tempDir, PURPLE_GERMLINE_GENE_DATA_CSV);
        String output = Files.readString(outputFile.toPath());
        List<String> lines = Arrays.asList(output.split("\n"));
        Assert.assertEquals(6, lines.size());
        checkLine(lines.get(1), "chr20", "1914001", "1915000", "DEL", "2");
        checkLine(lines.get(2), "chr20", "58692001", "58693000", "DEL", "3");
        checkLine(lines.get(3), "chr20", "58981001", "58982000", "DEL", "2");
    }

    private void checkLine(String line, String... expected)
    {
        List<String> parts = Arrays.asList(line.split(","));
        Assert.assertEquals(expected.length, parts.size());
        for(int i = 0; i < expected.length; i++)
        {
            Assert.assertEquals(expected[i], parts.get(i));
        }
    }
}
