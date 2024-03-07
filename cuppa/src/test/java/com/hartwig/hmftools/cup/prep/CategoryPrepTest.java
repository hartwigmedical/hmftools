package com.hartwig.hmftools.cup.prep;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.cup.feature.FeaturePrep;
import com.hartwig.hmftools.cup.rna.GeneExpressionPrep;
import com.hartwig.hmftools.cup.somatics.SomaticVariantPrep;
import com.hartwig.hmftools.cup.svs.StructuralVariantPrep;
import com.hartwig.hmftools.cup.traits.SampleTraitPrep;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class CategoryPrepTest
{
    public final String selectedSampleId = "COLO829v003T";

    public final PrepConfig prepConfig = new TestPrepConfigBuilder()
            .sampleIds(List.of(selectedSampleId))
            .categories(CategoryType.getAllCategories())
            .refGenomeVersion("V37")
            .sampleDataDir(TestPrepConfigBuilder.TEST_SAMPLE_DATA_DIR + "*")
            .outputDir("/tmp/")
            .build();

    @Test
    public void canExtractSnvFeatures()
    {
        SomaticVariantPrep prep = new SomaticVariantPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(selectedSampleId);

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
    }

    @Test
    public void canExtractSvFeatures()
    {
        StructuralVariantPrep prep = new StructuralVariantPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(selectedSampleId);

        assertEquals(6, dataItems.size());
        assertEquals(dataItems.get(0), new DataItem(DataSource.DNA, ItemType.SV_COUNT, "LINE", "3"));
        assertEquals(dataItems.get(1), new DataItem(DataSource.DNA, ItemType.SV_COUNT, "SIMPLE_DEL_20KB_1MB", "20"));
        assertEquals(dataItems.get(2), new DataItem(DataSource.DNA, ItemType.SV_COUNT, "SIMPLE_DUP_32B_200B", "3"));
        assertEquals(dataItems.get(3), new DataItem(DataSource.DNA, ItemType.SV_COUNT, "SIMPLE_DUP_100KB_5MB", "2"));
        assertEquals(dataItems.get(4), new DataItem(DataSource.DNA, ItemType.SV_COUNT, "MAX_COMPLEX_SIZE", "8"));
        assertEquals(dataItems.get(5), new DataItem(DataSource.DNA, ItemType.SV_COUNT, "TELOMERIC_SGL", "0"));
    }

    @Test
    public void canExtractTraitFeatures()
    {
        SampleTraitPrep prep = new SampleTraitPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(selectedSampleId);

        assertEquals(3, dataItems.size());
        assertEquals(dataItems.get(0), new DataItem(DataSource.DNA, ItemType.SAMPLE_TRAIT, "is_male", "1"));
        assertEquals(dataItems.get(1), new DataItem(DataSource.DNA, ItemType.TUMOR_MUTATIONAL_BURDEN, "indels_per_mb", "0.1207"));
        assertEquals(dataItems.get(2), new DataItem(DataSource.DNA, ItemType.SAMPLE_TRAIT, "whole_genome_duplication", "1"));
    }

    @Test
    public void canExtractEventFeatures()
    {
        FeaturePrep prep = new FeaturePrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(selectedSampleId);

        assertEquals(9, dataItems.size());
        assertEquals(dataItems.get(1), new DataItem(DataSource.DNA, ItemType.DRIVER, "BRAF.mut", "1.0"));
    }

    @Test
    public void canExtractGeneExpFeatures()
    {
        GeneExpressionPrep prep = new GeneExpressionPrep(prepConfig);
        List<DataItem> dataItems = prep.extractSampleData(selectedSampleId);

        assertEquals(2, dataItems.size());
        assertEquals(dataItems.get(0), new DataItem(DataSource.RNA, ItemType.EXPRESSION, "BRAF", "3.434e+00"));
        assertEquals(dataItems.get(1), new DataItem(DataSource.RNA, ItemType.EXPRESSION, "TP53", "3.870e+00"));
    }
}
