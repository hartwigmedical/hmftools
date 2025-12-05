package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.cup.TestPrepConfigBuilder.TEST_TUMOR_SAMPLE_ID;
import static com.hartwig.hmftools.cup.TestPrepConfigBuilder.TEST_TUMOR_SAMPLE_RNA_ID;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.cup.prep.CategoryType;
import com.hartwig.hmftools.cup.drivers.DriverPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.DataSource;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;
import com.hartwig.hmftools.cup.rna.AltSpliceJunctionPrep;
import com.hartwig.hmftools.cup.rna.GeneExpressionPrep;
import com.hartwig.hmftools.cup.somatics.SomaticVariantPrep;
import com.hartwig.hmftools.cup.svs.StructuralVariantPrep;
import com.hartwig.hmftools.cup.traits.SampleTraitPrep;

import org.junit.Test;

public class CategoryPrepTest
{
    public final PrepConfig prepConfig = new TestPrepConfigBuilder()
            .sampleIds(List.of(TEST_TUMOR_SAMPLE_ID))
            .rnaSampleIds(List.of(TEST_TUMOR_SAMPLE_RNA_ID))
            .categories(CategoryType.getAllCategories())
            .refGenomeVersion("V37")
            .sampleDataDir(TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*")
            .altSpliceJunctionSites(TestPrepConfigBuilder.TEST_ALT_SPLICE_JUNCTION_SITES)
            .outputDir("/tmp/")
            .build();

    private HashMap<String, String> makeDataItemsMap(List<DataItem> dataItems)
    {
        HashMap<String, String> dataItemsMap = new HashMap<>();
        for(DataItem dataItem : dataItems)
            dataItemsMap.put(dataItem.Index.Key, dataItem.Value);

        return dataItemsMap;
    }

    @Test
    public void canExtractSnvFeatures()
    {
        SomaticVariantPrep prep = new SomaticVariantPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(TEST_TUMOR_SAMPLE_ID);

        // Check that all required item types exist
        List<ItemType> itemTypesUnique = dataItems.stream().map(o -> o.Index.Type).distinct().collect(Collectors.toList());
        List<ItemType> itemTypesExpected = List.of(
                ItemType.SNV96,
                ItemType.TUMOR_MUTATIONAL_BURDEN,
                ItemType.GEN_POS,
                ItemType.SIGNATURE
        );

        for(ItemType itemTypeExpected : itemTypesExpected){
            assertTrue(itemTypesUnique.contains(itemTypeExpected));
        }

        // Check one value of each type
        HashMap<String, String> dataItemsMap = makeDataItemsMap(dataItems);
        assertEquals("2", dataItemsMap.get("C>T_TCC"));
        assertEquals("8", dataItemsMap.get("snv_count"));
        assertEquals("6.4", dataItemsMap.get("SIG_7_UV"));

        assertEquals("1", dataItemsMap.get("1_500000"));
        assertEquals("1", dataItemsMap.get("1_1000000"));
        assertEquals("1", dataItemsMap.get("7_55000000"));
    }

    @Test
    public void canExtractSnvFeaturesFromLiftoverGenericVariantsFiles()
    {
        PrepConfig prepConfig = new TestPrepConfigBuilder()
                .sampleIds(List.of(TEST_TUMOR_SAMPLE_ID))
                .categories(List.of(CategoryType.SNV))
                .refGenomeVersion("V38")
                .somaticVariantsDir(TestPrepConfigBuilder.TEST_SOMATIC_VARIANTS_DIR)
                .build();

        SomaticVariantPrep prep = new SomaticVariantPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(TEST_TUMOR_SAMPLE_ID);

        // Check one value of each type
        HashMap<String, String> dataItemsMap = makeDataItemsMap(dataItems);
        assertEquals("2", dataItemsMap.get("C>T_TCC"));
        assertEquals("8", dataItemsMap.get("snv_count"));
        assertEquals("6.4", dataItemsMap.get("SIG_7_UV"));

        assertEquals("1", dataItemsMap.get("chr1_500000"));
        assertEquals("1", dataItemsMap.get("chr1_1000000"));
        assertEquals("1", dataItemsMap.get("chr7_55000000"));
    }

    @Test
    public void canExtractSvFeatures()
    {
        StructuralVariantPrep prep = new StructuralVariantPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(TEST_TUMOR_SAMPLE_ID);

        assertEquals(6, dataItems.size());
        assertEquals(new DataItem(DataSource.DNA, ItemType.SV_COUNT, "LINE", "3"), dataItems.get(0));
        assertEquals(new DataItem(DataSource.DNA, ItemType.SV_COUNT, "SIMPLE_DEL_20KB_1MB", "25"), dataItems.get(1));
        assertEquals(new DataItem(DataSource.DNA, ItemType.SV_COUNT, "SIMPLE_DUP_32B_200B", "4"), dataItems.get(2));
        assertEquals(new DataItem(DataSource.DNA, ItemType.SV_COUNT, "SIMPLE_DUP_100KB_5MB", "2"), dataItems.get(3));
        assertEquals(new DataItem(DataSource.DNA, ItemType.SV_COUNT, "MAX_COMPLEX_SIZE", "8"), dataItems.get(4));
        assertEquals(new DataItem(DataSource.DNA, ItemType.SV_COUNT, "TELOMERIC_SGL", "0"), dataItems.get(5));
    }

    @Test
    public void canExtractTraitFeatures()
    {
        SampleTraitPrep prep = new SampleTraitPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(TEST_TUMOR_SAMPLE_ID);

        assertEquals(3, dataItems.size());
        assertEquals(new DataItem(DataSource.DNA, ItemType.SAMPLE_TRAIT, "is_male", "1"), dataItems.get(0));
        assertEquals(new DataItem(DataSource.DNA, ItemType.TUMOR_MUTATIONAL_BURDEN, "indels_per_mb", "0.1207"), dataItems.get(1));
        assertEquals(new DataItem(DataSource.DNA, ItemType.SAMPLE_TRAIT, "whole_genome_duplication", "1"), dataItems.get(2));
    }

    @Test
    public void canExtractEventFeatures()
    {
        DriverPrep prep = new DriverPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(TEST_TUMOR_SAMPLE_ID);

        assertEquals(8, dataItems.size());
        assertEquals(new DataItem(DataSource.DNA, ItemType.DRIVER, "BRAF.mut", "1.0000"), dataItems.get(0));
        assertEquals(new DataItem(DataSource.DNA, ItemType.DRIVER, "MYC.amp", "1.0000"), dataItems.get(1));
        assertEquals(new DataItem(DataSource.DNA, ItemType.DRIVER, "CDKN2A.mut", "1.0000"), dataItems.get(2));
        assertEquals(new DataItem(DataSource.DNA, ItemType.DRIVER, "SFTPB.indel", "1.0000"), dataItems.get(3));
        assertEquals(new DataItem(DataSource.DNA, ItemType.FUSION, "TMPRSS2_ERG", "1.0000"), dataItems.get(4));
        assertEquals(new DataItem(DataSource.DNA, ItemType.FUSION, "BRAF_PROM3", "1.0000"), dataItems.get(5));
        assertEquals(new DataItem(DataSource.DNA, ItemType.FUSION, "FGFR2_PROM5", "1.0000"), dataItems.get(6));
        assertEquals(new DataItem(DataSource.DNA, ItemType.VIRUS, "MERKEL", "1.0000"), dataItems.get(7));
    }

    @Test
    public void canExtractGeneExpFeatures()
    {
        GeneExpressionPrep prep = new GeneExpressionPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(TEST_TUMOR_SAMPLE_RNA_ID);

        assertEquals(2, dataItems.size());
        assertEquals(new DataItem(DataSource.RNA, ItemType.EXPRESSION, "BRAF", "3.434e+00"), dataItems.get(0));
        assertEquals(new DataItem(DataSource.RNA, ItemType.EXPRESSION, "TP53", "3.870e+00"), dataItems.get(1));
    }

    @Test
    public void canExtractAltSjFeatures()
    {
        AltSpliceJunctionPrep prep = new AltSpliceJunctionPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(TEST_TUMOR_SAMPLE_RNA_ID);

        assertEquals(4, dataItems.size());
        assertEquals(new DataItem(DataSource.RNA, ItemType.ALT_SJ, "7;140426316;140439612", "2"), dataItems.get(0));
        assertEquals(new DataItem(DataSource.RNA, ItemType.ALT_SJ, "10;89623836;89623905", "1"), dataItems.get(1));
        assertEquals(new DataItem(DataSource.RNA, ItemType.ALT_SJ, "13;32901736;32901958", "1"), dataItems.get(2));
        assertEquals(new DataItem(DataSource.RNA, ItemType.ALT_SJ, "20;46474688;46476232", "1"), dataItems.get(3));
    }
}
