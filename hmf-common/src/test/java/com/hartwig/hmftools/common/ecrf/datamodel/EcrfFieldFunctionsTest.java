package com.hartwig.hmftools.common.ecrf.datamodel;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class EcrfFieldFunctionsTest {
    @Test
    public void canConvertToName() {
        final String studyEventOID = "SE.Study";
        final String formOID = "FRM.Form";
        final String itemGroupOID = "GRP.ItemGroup";
        final String itemOID = "GRP.ItemGroup.Item";

        assertEquals("STUDY.FORM.ITEMGROUP.ITEM", EcrfFieldFunctions.name(studyEventOID, formOID, itemGroupOID, itemOID));
    }
}