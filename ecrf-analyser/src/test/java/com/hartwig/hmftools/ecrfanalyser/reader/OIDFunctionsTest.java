package com.hartwig.hmftools.ecrfanalyser.reader;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class OIDFunctionsTest {

    @Test
    public void canExtractEcrfFields() {
        String category = "cat";
        String fieldName = "fldName";

        String OID = OIDFunctions.toOID(category, fieldName);

        assertEquals(category, OIDFunctions.category(OID));
        assertEquals(fieldName, OIDFunctions.fieldName(OID));
    }

    @Test
    public void canExtractNoCategory() {
        String category = OIDFunctions.NO_CATEGORY;
        String fieldName = "fldName";

        String OID = OIDFunctions.toOID(category, fieldName);

        assertEquals(category, OIDFunctions.category(OID));
        assertEquals(fieldName, OIDFunctions.fieldName(OID));
    }
}