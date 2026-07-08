package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.generateFileName;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.read;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.toLines;
import static com.hartwig.hmftools.compar.common.FieldConfigFile.write;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
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

    private static String writeTempFile(final List<String> lines) throws IOException
    {
        File file = File.createTempFile("field_config", ".tsv");
        file.deleteOnExit();
        Files.write(file.toPath(), lines);
        return file.getPath();
    }

    @Test
    public void readParsesAllColumnsIntoFieldOverride() throws IOException
    {
        String filename = writeTempFile(List.of(HEADER, "PURITY\tDoubleField\tdouble\ttrue\t5.0\t20.0%"));

        List<FieldOverride> overrides = read(filename);

        assertEquals(1, overrides.size());
        FieldOverride override = overrides.get(0);
        assertEquals("PURITY", override.Category);
        assertEquals("DoubleField", override.Field);
        assertEquals("true", override.Compared);
        assertEquals("5.0", override.AbsoluteThreshold);
        assertEquals("20.0%", override.PercentThreshold);
    }

    @Test
    public void readPreservesBlankOptionalColumnsAsEmptyStrings() throws IOException
    {
        String filename = writeTempFile(List.of(HEADER, "PURITY\tDoubleField\tdouble\t\t\t"));

        FieldOverride override = read(filename).get(0);

        assertEquals("", override.Compared);
        assertEquals("", override.AbsoluteThreshold);
        assertEquals("", override.PercentThreshold);
    }

    @Test
    public void readHandlesMultipleRowsInFileOrder() throws IOException
    {
        String filename = writeTempFile(List.of(
                HEADER,
                "PURITY\tPurityField\tstring\ttrue\tnone\tnone",
                "DRIVER\tDriverField\tstring\tfalse\tnone\tnone"));

        List<FieldOverride> overrides = read(filename);

        assertEquals(2, overrides.size());
        assertEquals("PurityField", overrides.get(0).Field);
        assertEquals("DriverField", overrides.get(1).Field);
    }

    @Test
    public void readRoundTripsFileWrittenByWrite() throws IOException
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new DoubleField("DoubleField", i -> 0d, false, 5.0, 0.25, "%.1f"));

        File file = File.createTempFile("field_config", ".tsv");
        file.deleteOnExit();

        write(file.getPath(), fieldConfig, Set.of(PURITY));

        FieldOverride override = read(file.getPath()).get(0);

        assertEquals("PURITY", override.Category);
        assertEquals("DoubleField", override.Field);
        assertEquals("false", override.Compared);
        assertEquals("5.0", override.AbsoluteThreshold);
        assertEquals("25.0%", override.PercentThreshold);
    }
}
