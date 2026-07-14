package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.compar.common.field.BooleanField;
import com.hartwig.hmftools.compar.common.field.DisplayOnlyField;
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

        FieldOverride override = new FieldOverride(PURITY.toString(), "BoolField", "false", "", "");
        fieldConfig.applyOverrides(List.of(override), false);

        assertFalse(getField(fieldConfig, PURITY, "BoolField").isCompared());
    }

    @Test
    public void applyOverrideUpdatesThresholds()
    {
        FieldConfig fieldConfig = new FieldConfig();
        DoubleField beforeOverride = new DoubleField("DoubleField", i -> 0d, true, 0.1,
                null, "%.2f");
        fieldConfig.registerField(PURITY, beforeOverride);

        FieldOverride override = new FieldOverride(PURITY.toString(), "DoubleField", "", "0.5", "20%");
        fieldConfig.applyOverrides(List.of(override), false);

        Field afterOverride = getField(fieldConfig, PURITY, "DoubleField");
        assertTrue(afterOverride.isCompared());
        assertEquals(0.5, afterOverride.absoluteThreshold(), 1e-9);
        assertEquals(0.2, afterOverride.percentThreshold(), 1e-9);
    }

    @Test
    public void applyOverrideAcceptsFractionForPercentThreshold()
    {
        FieldConfig fieldConfig = new FieldConfig();
        DoubleField beforeOverride = new DoubleField("DoubleField", i -> 0d, true, null,
                null, "%.2f");
        fieldConfig.registerField(PURITY, beforeOverride);

        FieldOverride override = new FieldOverride(PURITY.toString(), "DoubleField", "", "", "0.5");
        fieldConfig.applyOverrides(List.of(override), false);

        assertEquals(0.5, getField(fieldConfig, PURITY, "DoubleField").percentThreshold(), 1e-9);
    }

    @Test
    public void applyOverrideNoneClearsThresholds()
    {
        FieldConfig fieldConfig = new FieldConfig();
        DoubleField beforeOverride = new DoubleField("DoubleField", i -> 0d, true, 0.1,
                0.2, "%.2f");
        fieldConfig.registerField(PURITY, beforeOverride);

        FieldOverride override = new FieldOverride(PURITY.toString(), "DoubleField", "", "none", "none");
        fieldConfig.applyOverrides(List.of(override), false);

        Field afterOverride = getField(fieldConfig, PURITY, "DoubleField");
        assertNull(afterOverride.absoluteThreshold());
        assertNull(afterOverride.percentThreshold());
    }

    @Test
    public void applyOverrideLeavesUnspecifiedSettingsUnchanged()
    {
        FieldConfig fieldConfig = new FieldConfig();
        DoubleField beforeOverride = new DoubleField("DoubleField", i -> 0d, true, 0.1,
                0.2, "%.2f");
        fieldConfig.registerField(PURITY, beforeOverride);

        fieldConfig.applyOverrides(List.of(new FieldOverride(PURITY.toString(), "DoubleField", "", "", "")), false);

        Field afterOverride = getField(fieldConfig, PURITY, "DoubleField");
        assertTrue(afterOverride.isCompared());
        assertEquals(0.1, afterOverride.absoluteThreshold(), 1e-9);
        assertEquals(0.2, afterOverride.percentThreshold(), 1e-9);
    }

    @Test
    public void applyOverrideIgnoresUnknownCategory()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new BooleanField("BoolField", i -> true, true));

        FieldOverride override = new FieldOverride("NOT_A_CATEGORY", "BoolField", "false", "", "");
        fieldConfig.applyOverrides(List.of(override), false);

        assertTrue(getField(fieldConfig, PURITY, "BoolField").isCompared());
    }

    @Test
    public void applyOverrideIgnoresUnknownFieldName()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new BooleanField("BoolField", i -> true, true));

        FieldOverride override = new FieldOverride(PURITY.toString(), "NotARealField", "false", "", "");
        fieldConfig.applyOverrides(List.of(override), false);

        assertTrue(getField(fieldConfig, PURITY, "BoolField").isCompared());
    }

    @Test
    public void applyOverrideOnUnsupportedSettingIsNoOpForThatSetting()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new BooleanField("BoolField", i -> true, true));

        // boolean fields don't support thresholds - should record an error and leave the field otherwise intact
        FieldOverride override = new FieldOverride(PURITY.toString(), "BoolField", "false", "5.0", "");
        fieldConfig.applyOverrides(List.of(override), false);

        Field afterOverride = getField(fieldConfig, PURITY, "BoolField");
        assertFalse(afterOverride.isCompared());
        assertNull(afterOverride.absoluteThreshold());
        assertEquals(1, fieldConfig.errorMessages().size());
    }

    @Test
    public void applyOverrideOnlyAffectsTargetedCategoryAndField()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new BooleanField("BoolField", i -> true, true));
        fieldConfig.registerField(DRIVER, new BooleanField("BoolField", i -> true, true));

        FieldOverride override = new FieldOverride(PURITY.toString(), "BoolField", "false", "", "");
        fieldConfig.applyOverrides(List.of(override), false);

        assertFalse(getField(fieldConfig, PURITY, "BoolField").isCompared());
        assertTrue(getField(fieldConfig, DRIVER, "BoolField").isCompared());
    }

    @Test
    public void applyOverridesStrictModeRecordsErrorForMissingColumns()
    {
        FieldConfig fieldConfig = new FieldConfig();
        DoubleField beforeOverride = new DoubleField("DoubleField", i -> 0d, true, 0.1,
                0.2, "%.2f");
        fieldConfig.registerField(PURITY, beforeOverride);

        FieldOverride override = new FieldOverride(PURITY.toString(), "DoubleField", "", "", "");
        fieldConfig.applyOverrides(List.of(override), true);

        assertEquals(3, fieldConfig.errorMessages().size());
        assertTrue(fieldConfig.warnings().isEmpty());

        Field field = getField(fieldConfig, PURITY, "DoubleField");
        assertTrue(field.isCompared());
        assertEquals(0.1, field.absoluteThreshold(), 1e-9);
        assertEquals(0.2, field.percentThreshold(), 1e-9);
    }

    @Test
    public void applyOverridesRecordsErrorForMalformedThreshold()
    {
        for(boolean strictFieldConfig : List.of(false, true))
        {
            FieldConfig fieldConfig = new FieldConfig();
            DoubleField beforeOverride = new DoubleField("DoubleField", i -> 0d, true, 0.1,
                    0.2, "%.2f");
            fieldConfig.registerField(PURITY, beforeOverride);

            FieldOverride override = new FieldOverride(PURITY.toString(), "DoubleField", "true",
                    "not-a-number", "10%");
            fieldConfig.applyOverrides(List.of(override), strictFieldConfig);

            assertEquals(1, fieldConfig.errorMessages().size());

            Field field = getField(fieldConfig, PURITY, "DoubleField");
            assertNull(field.absoluteThreshold());
        }
    }

    @Test
    public void applyOverridesRecordsErrorForMalformedCompared()
    {
        for(boolean strictFieldConfig : List.of(false, true))
        {
            FieldConfig fieldConfig = new FieldConfig();
            fieldConfig.registerField(PURITY, new BooleanField("BoolField", i -> true, true));

            FieldOverride override = new FieldOverride(PURITY.toString(), "BoolField", "yes", "none",
                    "none");
            fieldConfig.applyOverrides(List.of(override), strictFieldConfig);

            assertEquals(1, fieldConfig.errorMessages().size());
            assertTrue(getField(fieldConfig, PURITY, "BoolField").isCompared());
        }
    }

    @Test
    public void applyOverridesStrictModeFlagsMissingOverride()
    {
        FieldConfig fieldConfig = new FieldConfig();
        fieldConfig.registerField(PURITY, new BooleanField("OverriddenField", i -> true, true));
        fieldConfig.registerField(PURITY, new BooleanField("MissingField", i -> true, true));
        fieldConfig.registerField(PURITY, new DisplayOnlyField("DisplayOnlyField", i -> "", i -> true));

        FieldOverride override = new FieldOverride(PURITY.toString(), "OverriddenField", "false", "none",
                "none");
        fieldConfig.applyOverrides(List.of(override), true);

        assertEquals(1, fieldConfig.errorMessages().size());
        assertTrue(fieldConfig.errorMessages().get(0).contains("MissingField"));
    }
}
