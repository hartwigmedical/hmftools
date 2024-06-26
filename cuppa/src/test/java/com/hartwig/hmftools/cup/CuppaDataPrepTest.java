package com.hartwig.hmftools.cup;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.cup.prep.CategoryType;
import com.hartwig.hmftools.cup.prep.CuppaDataPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.DataItemMatrix;
import com.hartwig.hmftools.cup.prep.DataItemsIO;
import com.hartwig.hmftools.cup.prep.DataSource;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;

import org.apache.commons.io.FileUtils;

import org.junit.Test;

public class CuppaDataPrepTest
{
    File TMP_DIR = new File(System.getProperty("java.io.tmpdir") + "/CuppaDataPrepTest/");

    @Test
    public void canRunSingleSamplePrep() throws IOException
    {
        TMP_DIR.mkdir();

        String selectedSampleId = "TUMOR_SAMPLE";

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

        assertEquals(6228,  dataItems.size());

        assertEquals(new DataItem(DataSource.DNA, ItemType.SNV96, "C>T_TCC", "2"), dataItems.get(45));
        assertEquals(new DataItem(DataSource.RNA, ItemType.ALT_SJ, "13;32901736;32901958", "1"), dataItems.get(dataItems.size()-1));

        // Check output file exists
        File outputPath = new File(cuppaDataPrep.getOutputPath(null));
        assertTrue(outputPath.exists());

        FileUtils.deleteDirectory(TMP_DIR);
    }

    @Test
    public void canRunMultiSamplePrep() throws IOException
    {
        TMP_DIR.mkdir();

        PrepConfig prepConfig = new TestPrepConfigBuilder()
                .sampleIds(Arrays.asList(
                        "MINIMAL_SAMPLE", "MINIMAL_SAMPLE", "MINIMAL_SAMPLE",
                        "MINIMAL_SAMPLE_NO_RNA",
                        "TUMOR_SAMPLE"
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

        HashMap<CategoryType, DataItemMatrix> dataItemMatricesByCategory = cuppaDataPrep.mDataItemMatricesByCategory;
        for(CategoryType categoryType : categoryTypes)
        {
            // Check that values are exactly the same between duplicate "MINIMAL_SAMPLE" samples.
            // Thread unsafe operations lead to different feature values even though the samples are the same.
            DataItemMatrix dataItemMatrix = dataItemMatricesByCategory.get(categoryType);

            assertEquals(dataItemMatrix.getSampleFeatureValues(0), dataItemMatrix.getSampleFeatureValues(1));

            assertEquals(dataItemMatrix.getSampleFeatureValues(2), dataItemMatrix.getSampleFeatureValues(3));

            // Check that values of the "MINIMAL_SAMPLE" and "TUMOR_SAMPLE" are different
            assertNotEquals(dataItemMatrix.getSampleFeatureValues(0), dataItemMatrix.getSampleFeatureValues(4));

            // Check output files exist
            File outputFile = new File(cuppaDataPrep.getOutputPath(categoryType));
            assertTrue(outputFile.exists());
        }

        FileUtils.deleteDirectory(TMP_DIR);
    }

    @Test
    public void canReadSingleSampleOutputFile()
    {
        String path = Resources.getResource("prep_output/TUMOR_SAMPLE.cuppa_data.modified.tsv").getPath();
        List<DataItem> dataItems = DataItemsIO.readDataItemList(path);

        assertEquals(new DataItem(DataSource.DNA, ItemType.SNV96, "C>T_ACG", "1"), dataItems.get(0));
        assertEquals(new DataItem(DataSource.DNA, ItemType.TUMOR_MUTATIONAL_BURDEN, "snv_count", "8"), dataItems.get(1));
        assertEquals(new DataItem(DataSource.DNA, ItemType.SIGNATURE, "SIG_7_UV", "6.4"), dataItems.get(2));
        assertEquals(new DataItem(DataSource.DNA, ItemType.GEN_POS, "1_500000", "1"), dataItems.get(3));
        assertEquals(new DataItem(DataSource.DNA, ItemType.SV_COUNT, "LINE", "3"), dataItems.get(4));
        assertEquals(new DataItem(DataSource.DNA, ItemType.SAMPLE_TRAIT, "is_male", "1"), dataItems.get(5));
        assertEquals(new DataItem(DataSource.DNA, ItemType.DRIVER, "BRAF.mut", "1"), dataItems.get(6));
        assertEquals(new DataItem(DataSource.DNA, ItemType.FUSION, "TMPRSS2_ERG", "1"), dataItems.get(7));
        assertEquals(new DataItem(DataSource.DNA, ItemType.VIRUS, "MERKEL", "1"), dataItems.get(8));
        assertEquals(new DataItem(DataSource.RNA, ItemType.EXPRESSION, "BRAF", "3.43E+00"), dataItems.get(9));
        assertEquals(new DataItem(DataSource.RNA, ItemType.ALT_SJ, "7;140426316;140439612", "2"), dataItems.get(10));
    }

    @Test
    public void canReadMultiSampleOutputFile()
    {
        String path = Resources.getResource("prep_output/cuppa_data.cohort.gene_exp.modified.tsv").getPath();
        DataItemMatrix dataItemMatrix = DataItemsIO.readDataItemMatrix(path);

        List<DataItem.Index> matrixIndexes =  dataItemMatrix.Indexes;

        assertEquals(new DataItem.Index(DataSource.RNA, ItemType.EXPRESSION, "BRAF"), matrixIndexes.get(0));
        assertEquals(new String[] { null, "3.43E+00" }, dataItemMatrix.get(matrixIndexes.get(0)));

        assertEquals(new DataItem.Index(DataSource.RNA, ItemType.EXPRESSION, "TP53"), matrixIndexes.get(1));
        assertEquals(new String[] { null, "3.87E+00" }, dataItemMatrix.get(matrixIndexes.get(1)));
    }
}
