package com.hartwig.hmftools.ecrfanalyser.reader;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ItemDefToEcrfFieldTest {

    @Test
    public void canExtractEcrfFields() {
        String category = "cat";
        String fieldName = "fldName";
        String description = "desc";
        ItemDef def = new ItemDef(ItemDefTestFunctions.toOID(category, fieldName), description, null);

        assertEquals(category, ItemDefToEcrfField.category(def));
        assertEquals(fieldName, ItemDefToEcrfField.fieldName(def));
        assertEquals(description, ItemDefToEcrfField.description(def));
    }
}