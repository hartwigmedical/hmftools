package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class OIDFunctionsTest {

    @Test
    public void canConvertBetweenFieldNameAndOID() {
        String ecrfFieldName = "cat.fldName";
        String OID = OIDFunctions.toOID(ecrfFieldName);

        assertEquals(ecrfFieldName, OIDFunctions.toEcrfFieldName(OID));
    }

    @Test
    public void canExtractNoCategory() {
        String name = OIDFunctions.NO_CATEGORY + OIDFunctions.OID_SEPARATOR + "fldName";
        String OID = OIDFunctions.toOID(name);

        assertEquals(name, OIDFunctions.toEcrfFieldName(OID));
    }
}