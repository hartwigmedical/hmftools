package com.hartwig.hmftools.cup.prep;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.cuppa.CategoryType;

import org.junit.Test;

public class CuppaDataPrepTest
{
    File TMP_DIR = new File(System.getProperty("java.io.tmpdir") + "/CuppaDataPrepTest/");

    private void deleteTmpDir()
    {
        for(File file : TMP_DIR.listFiles())
            file.delete();

        TMP_DIR.delete();
    }

    @Test
    public void canRunSingleSamplePrep()
    {
        TMP_DIR.mkdir();

        String selectedSampleId = "COLO829v003T";

        PrepConfig prepConfig = new TestPrepConfigBuilder()
                .sampleIds(List.of(selectedSampleId))
                .categories(CategoryType.getDnaCategories())
                .refGenomeVersion("V37")
                .sampleDataDir(TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*")
                .outputDir(TMP_DIR.toString())
                .build();

        assertTrue(prepConfig.isSingleSample());

        CuppaDataPrep cuppaDataPrep = new CuppaDataPrep(prepConfig);
        cuppaDataPrep.run();

        File outputPath = new File(cuppaDataPrep.getOutputPath(null));
        assertTrue(outputPath.exists());

        deleteTmpDir();
    }

    @Test
    public void canRunMultiSamplePrep()
    {
        TMP_DIR.mkdir();

        PrepConfig prepConfig = new TestPrepConfigBuilder()
                .sampleIds(Arrays.asList("COLO829v003T", "COLO829v003T_modified"))
                .categories(CategoryType.getDnaCategories())
                .refGenomeVersion("V37")
                .sampleDataDir(TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*")
                .outputDir(TMP_DIR.toString())
                .threads(4)
                .build();

        CuppaDataPrep cuppaDataPrep = new CuppaDataPrep(prepConfig);
        cuppaDataPrep.run();

        for(CategoryType categoryType : CategoryType.getDnaCategories())
        {
            File outputFile = new File(cuppaDataPrep.getOutputPath(categoryType));
            assertTrue(outputFile.exists());
        }

        deleteTmpDir();
    }
}
