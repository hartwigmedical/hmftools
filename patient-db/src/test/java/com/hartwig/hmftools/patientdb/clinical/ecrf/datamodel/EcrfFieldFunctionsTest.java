package com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class EcrfFieldFunctionsTest {

    @Test
    public void canConvertToName() {
        String studyEventOID = "SE.Study";
        String formOID = "FRM.Form";
        String itemGroupOID = "GRP.ItemGroup";
        String itemOID = "GRP.ItemGroup.Item";

        assertEquals("STUDY.FORM.ITEMGROUP.ITEM", EcrfFieldFunctions.name(studyEventOID, formOID, itemGroupOID, itemOID));
    }
}