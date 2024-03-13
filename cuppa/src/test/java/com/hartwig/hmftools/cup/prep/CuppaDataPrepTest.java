package com.hartwig.hmftools.cup.prep;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
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
                .categories(CategoryType.getAllCategories())
                .refGenomeVersion("V37")
                .sampleDataDir(TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*")
                .outputDir(TMP_DIR.toString())
                .build();

        assertTrue(prepConfig.isSingleSample());

        CuppaDataPrep cuppaDataPrep = new CuppaDataPrep(prepConfig);
        cuppaDataPrep.run(true);

        // Check values
        List<DataItem> dataItems = cuppaDataPrep.mDataItems;

        assertEquals(6348,  dataItems.size());

        assertEquals(
                new DataItem(DataSource.DNA, ItemType.SNV96, "C>A_ACA", "133"),
                dataItems.get(0)
        );

        assertEquals(
                new DataItem(DataSource.RNA, ItemType.ALT_SJ, "13;32901736;32901958", "1"),
                dataItems.get(dataItems.size()-1)
        );

        // Check output file exists
        File outputPath = new File(cuppaDataPrep.getOutputPath(null));
        assertTrue(outputPath.exists());
        deleteTmpDir();
    }

    @Test
    public void canRunMultiSamplePrep()
    {
        TMP_DIR.mkdir();

        PrepConfig prepConfig = new TestPrepConfigBuilder()
                .sampleIds(Arrays.asList(
                        "COLO829v003T", "COLO829v003T",
                        "COLO829v003T_modified", "COLO829v003T_modified", "COLO829v003T_modified", "COLO829v003T_modified"
                ))
                .categories(CategoryType.getAllCategories())
                .refGenomeVersion("V37")
                .sampleDataDir(TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*")
                .outputDir(TMP_DIR.toString())
                .threads(5)
                .build();

        CuppaDataPrep cuppaDataPrep = new CuppaDataPrep(prepConfig);
        cuppaDataPrep.run(true);

        List<CategoryType> categoryTypes = cuppaDataPrep.mConfig.Categories;

        HashMap<CategoryType, DataItemMatrix>  dataItemMatricesByCategory = cuppaDataPrep.mDataItemMatricesByCategory;
        for(CategoryType categoryType : categoryTypes)
        {
            // Check that values of the 2 COLO829v003T runs are exactly the same.
            // Thread unsafe operations lead to different feature values even though the samples are the same.
            DataItemMatrix dataItemMatrix = dataItemMatricesByCategory.get(categoryType);
            assertEquals(
                    dataItemMatrix.getFeatureValuesBySampleIndex(0),
                    dataItemMatrix.getFeatureValuesBySampleIndex(1)
            );

            assertEquals(
                    dataItemMatrix.getFeatureValuesBySampleIndex(2),
                    dataItemMatrix.getFeatureValuesBySampleIndex(3)
            );

            // Check that values of the COLO829v003T and COLO829v003T_modified are different
            assertNotEquals(
                    dataItemMatrix.getFeatureValuesBySampleIndex(0),
                    dataItemMatrix.getFeatureValuesBySampleIndex(2)
            );

            // Check output files exist
            File outputFile = new File(cuppaDataPrep.getOutputPath(categoryType));
            assertTrue(outputFile.exists());
        }

        deleteTmpDir();
    }
}
