package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.compar.common.field.BooleanField;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;

import org.junit.Test;

public class FieldConfigTest
{
    private static Field getField(final FieldConfig fieldConfig, final CategoryType category, final String name)
    {
        return fieldConfig.getFields(category, List.of(name)).get(0);
    }

    @Test
    public void applyOverrideUpdatesCompared()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new BooleanField("BoolField", i -> true, true));

        fieldConfig.applyOverride(new FieldOverride(PURITY.toString(), "BoolField", "false", "", ""));

        assertFalse(getField(fieldConfig, PURITY, "BoolField").isCompared());
    }

    @Test
    public void applyOverrideUpdatesThresholds()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new DoubleField("DoubleField", i -> 0d, true, 0.1, null, "%.2f"));

        fieldConfig.applyOverride(new FieldOverride(PURITY.toString(), "DoubleField", "", "0.5", "20%"));

        Field field = getField(fieldConfig, PURITY, "DoubleField");
        assertTrue(field.isCompared());
        assertEquals(0.5, field.absoluteThreshold(), 1e-9);
        assertEquals(0.2, field.percentThreshold(), 1e-9);
    }

    @Test
    public void applyOverrideAcceptsFractionForPercentThreshold()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new DoubleField("DoubleField", i -> 0d, true, null, null, "%.2f"));

        fieldConfig.applyOverride(new FieldOverride(PURITY.toString(), "DoubleField", "", "", "0.5"));

        assertEquals(0.5, getField(fieldConfig, PURITY, "DoubleField").percentThreshold(), 1e-9);
    }

    @Test
    public void applyOverrideNoneClearsThresholds()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new DoubleField("DoubleField", i -> 0d, true, 0.1, 0.2, "%.2f"));

        fieldConfig.applyOverride(new FieldOverride(PURITY.toString(), "DoubleField", "", "none", "none"));

        Field field = getField(fieldConfig, PURITY, "DoubleField");
        assertNull(field.absoluteThreshold());
        assertNull(field.percentThreshold());
    }

    @Test
    public void applyOverrideLeavesUnspecifiedSettingsUnchanged()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new DoubleField("DoubleField", i -> 0d, true, 0.1, 0.2, "%.2f"));

        fieldConfig.applyOverride(new FieldOverride(PURITY.toString(), "DoubleField", "", "", ""));

        Field field = getField(fieldConfig, PURITY, "DoubleField");
        assertTrue(field.isCompared());
        assertEquals(0.1, field.absoluteThreshold(), 1e-9);
        assertEquals(0.2, field.percentThreshold(), 1e-9);
    }

    @Test
    public void applyOverrideIgnoresUnknownCategory()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new BooleanField("BoolField", i -> true, true));

        fieldConfig.applyOverride(new FieldOverride("NOT_A_CATEGORY", "BoolField", "false", "", ""));

        assertTrue(getField(fieldConfig, PURITY, "BoolField").isCompared());
    }

    @Test
    public void applyOverrideIgnoresUnknownFieldName()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new BooleanField("BoolField", i -> true, true));

        fieldConfig.applyOverride(new FieldOverride(PURITY.toString(), "NotARealField", "false", "", ""));

        assertTrue(getField(fieldConfig, PURITY, "BoolField").isCompared());
    }

    @Test
    public void applyOverrideOnUnsupportedSettingIsNoOpForThatSetting()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new BooleanField("BoolField", i -> true, true));

        // boolean fields don't support thresholds - should warn and leave the field otherwise intact
        fieldConfig.applyOverride(new FieldOverride(PURITY.toString(), "BoolField", "false", "5.0", ""));

        Field field = getField(fieldConfig, PURITY, "BoolField");
        assertFalse(field.isCompared());
        assertNull(field.absoluteThreshold());
    }

    @Test
    public void applyOverrideOnlyAffectsTargetedCategoryAndField()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new BooleanField("BoolField", i -> true, true));
        fieldConfig.registerField(DRIVER, new BooleanField("BoolField", i -> true, true));

        fieldConfig.applyOverride(new FieldOverride(PURITY.toString(), "BoolField", "false", "", ""));

        assertFalse(getField(fieldConfig, PURITY, "BoolField").isCompared());
        assertTrue(getField(fieldConfig, DRIVER, "BoolField").isCompared());
    }
}