package com.hartwig.hmftools.cup;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import com.hartwig.hmftools.cup.prep.CategoryType;
import com.hartwig.hmftools.cup.prep.PrepConfig;

import org.junit.Test;

public class PrepConfigTest
{
    private final String selectedSampleId = "sample_1";

    @Test
    public void canParseCommandLineArgsMultiSample()
    {
        String[] args = {
                // "-sample", sample,
                "-sample_id_file", TestPrepConfigBuilder.TEST_SAMPLE_ID_FILE,
                "-categories", "DNA",
                "-ref_genome_version", "V37",

                // Wild-cards are required to get subdirs correctly
                "-sample_data_dir", TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*",
                "-purple_dir", TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*/purple/",
                "-linx_dir", TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*/linx/",
                "-virus_dir", TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*/virus_interpreter/",

                "-output_dir", "/tmp/",
                "-isofox_dir", "/tmp/",
                "-ref_alt_sj_sites", TestPrepConfigBuilder.TEST_ALT_SPLICE_JUNCTION_SITES,
                "-threads", TestPrepConfigBuilder.TEST_THREADS.toString(),

                "-write_by_category"
        };

        PrepConfig prepConfig = TestPrepConfigBuilder.fromArgs(args);

        assertEquals(2, prepConfig.SampleIds.size());
        assertEquals("V37", prepConfig.RefGenVersion.toString());
        assertTrue(prepConfig.WriteByCategory);
        assertEquals((int) TestPrepConfigBuilder.TEST_THREADS, prepConfig.Threads);

        assertEquals(TestPrepConfigBuilder.TEST_ALT_SPLICE_JUNCTION_SITES, prepConfig.AltSpliceJunctionSites);

        String expectedPurpleDir = TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + selectedSampleId + "/purple/";
        String actualPurpleDir = prepConfig.getPurpleDataDir(selectedSampleId);
        assertEquals(expectedPurpleDir, actualPurpleDir);

        String expectedLinxDir = TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + selectedSampleId + "/linx/";
        String actualLinxDir = prepConfig.getLinxDataDir(selectedSampleId);
        assertEquals(expectedLinxDir, actualLinxDir);

        String expectedVirusDir = TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + selectedSampleId + "/virus_interpreter/";
        String actualVirusDir = prepConfig.getVirusDataDir(selectedSampleId);
        assertEquals(expectedVirusDir, actualVirusDir);
    }

    @Test
    public void sampleDirUsedWhenToolDirNotSpecified()
    {
        String[] args = {
                "-sample_id_file", TestPrepConfigBuilder.TEST_SAMPLE_ID_FILE,
                "-sample_data_dir", TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*",
        };

        PrepConfig prepConfig = TestPrepConfigBuilder.fromArgs(args);

        String expectedPurpleDir = TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + selectedSampleId;
        String actualPurpleDir = prepConfig.getPurpleDataDir(selectedSampleId);
        assertEquals(expectedPurpleDir, actualPurpleDir);
    }

    @Test
    public void canBuildFromTestPrepConfigBuilder()
    {
        PrepConfig prepConfig = new TestPrepConfigBuilder()
                .sampleIds(Arrays.asList("sample_1", "sample_2"))
                .categories(CategoryType.getDnaCategories())
                .refGenomeVersion("V37")
                .sampleDataDir(TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*")
                .linxDir(TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*/linx/")
                .outputDir("/Users/lnguyen/Desktop/test_output/")
                .build();

        String expectedPurpleDir = TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + selectedSampleId;
        String actualPurpleDir = prepConfig.getPurpleDataDir(selectedSampleId);
        assertEquals(expectedPurpleDir, actualPurpleDir);

        String expectedLinxDir = TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + selectedSampleId + "/linx/";
        String actualLinxDir = prepConfig.getLinxDataDir(selectedSampleId);
        assertEquals(expectedLinxDir, actualLinxDir);
    }
}
