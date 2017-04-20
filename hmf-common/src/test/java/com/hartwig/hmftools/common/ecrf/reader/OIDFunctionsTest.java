package com.hartwig.hmftools.common.ecrf.reader;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class OIDFunctionsTest {

    @Test
    public void canConvertBetweenFieldNameAndOID() {
        final String ecrfFieldName = "cat.fldName";
        final String OID = OIDFunctions.toOID(ecrfFieldName);

        assertEquals(ecrfFieldName, OIDFunctions.toEcrfFieldName(OID));
    }

    @Test
    public void canExtractNoCategory() {
        final String name = OIDFunctions.NO_CATEGORY + OIDFunctions.OID_SEPARATOR + "fldName";
        final String OID = OIDFunctions.toOID(name);

        assertEquals(name, OIDFunctions.toEcrfFieldName(OID));
    }
}