package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.generateFileName;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.toLines;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.compar.common.field.DisplayField;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.StringField;

import org.junit.Test;

public class FieldConfigFileTest
{
    private static final String HEADER = "category\tfield\tfieldType\tcompared\tabsoluteThreshold\tpercentThreshold";

    @Test
    public void generateFileNameAppendsFileNameWithDirSeparator()
    {
        assertEquals("dir" + File.separator + "field.config.compar.tsv", generateFileName("dir"));
        assertEquals("dir" + File.separator + "field.config.compar.tsv", generateFileName("dir" + File.separator));
    }

    @Test
    public void toLinesReturnsOnlyHeaderWhenNoCategoriesRequested()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new StringField("Name", i -> "", true));

        List<String> lines = toLines(fieldConfig, Set.of());

        assertEquals(List.of(HEADER), lines);
    }

    @Test
    public void toLinesOnlyIncludesRequestedCategories()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new StringField("PurityField", i -> "", true));
        fieldConfig.registerField(DRIVER, new StringField("DriverField", i -> "", true));

        List<String> lines = toLines(fieldConfig, Set.of(PURITY));

        assertEquals(List.of(HEADER, "PURITY\tPurityField\tstring\ttrue\tnone\tnone"), lines);
    }

    @Test
    public void toLinesExcludesDisplayFields()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new StringField("PurityField", i -> "", true));
        fieldConfig.registerField(PURITY, new DisplayField("DisplayOnly", i -> "", i -> true));

        List<String> lines = toLines(fieldConfig, Set.of(PURITY));

        assertEquals(List.of(HEADER, "PURITY\tPurityField\tstring\ttrue\tnone\tnone"), lines);
    }

    @Test
    public void toLinesOrdersFieldsByCategoryTypeEnumOrderRegardlessOfRegistrationOrder()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(DRIVER, new StringField("DriverField", i -> "", true));
        fieldConfig.registerField(PURITY, new StringField("PurityField", i -> "", true));

        List<String> lines = toLines(fieldConfig, Set.of(DRIVER, PURITY));

        assertEquals(List.of(
                HEADER,
                "PURITY\tPurityField\tstring\ttrue\tnone\tnone",
                "DRIVER\tDriverField\tstring\ttrue\tnone\tnone"), lines);
    }

    @Test
    public void toLinesFormatsIsComparedAndThresholds()
    {
        double absoluteThreshold = 5.0;
        double percentThreshold = 0.25;

        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY,
                new DoubleField("DoubleField", i -> 0d, false, absoluteThreshold, percentThreshold, "%.1f"));

        List<String> lines = toLines(fieldConfig, Set.of(PURITY));

        String expectedLine = "PURITY\tDoubleField\tdouble\tfalse\t" + absoluteThreshold + "\t" + (percentThreshold * 100) + "%";
        assertEquals(List.of(HEADER, expectedLine), lines);
    }

    @Test
    public void toLinesUsesNoneForNullThresholds()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new DoubleField("DoubleField", i -> 0d, true, null, null, "%.1f"));

        List<String> lines = toLines(fieldConfig, Set.of(PURITY));

        assertEquals(List.of(HEADER, "PURITY\tDoubleField\tdouble\ttrue\tnone\tnone"), lines);
    }
}