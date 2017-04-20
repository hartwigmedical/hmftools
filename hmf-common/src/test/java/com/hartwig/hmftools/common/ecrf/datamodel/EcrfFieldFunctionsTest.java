package com.hartwig.hmftools.common.ecrf.datamodel;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class EcrfFieldFunctionsTest {

    @Test
    public void canDetermineFieldRelevance() {
        verifyItemGroup(true, "Relevant");
        verifyItemGroup(false, "MetaData");
        verifyItem(false, "GROUP");
        verifyItem(false, "group3");
        verifyItem(true, "Relevant");
    }

    private static void verifyItemGroup(boolean isRelevant, @NotNull String itemGroupOID) {
        final EcrfField field = new EcrfField("SE.Study", "FRM.Form", "GRP." + itemGroupOID, "GRP.Item", "",
                Maps.newHashMap());
        assertEquals(isRelevant, EcrfFieldFunctions.isRelevant(field));
    }

    private static void verifyItem(boolean expected, @NotNull String itemOID) {
        final EcrfField field = new EcrfField("SE.Study", "FRM.Form", "GRP.ItemGroup", "GRP." + itemOID, "",
                Maps.newHashMap());
        assertEquals(expected, EcrfFieldFunctions.isRelevant(field));
    }

    @Test
    public void canConvertToName() {
        final String studyEventOID = "SE.Study";
        final String formOID = "FRM.Form";
        final String itemGroupOID = "GRP.ItemGroup";
        final String itemOID = "GRP.ItemGroup.Item";

        assertEquals("STUDY.FORM.ITEMGROUP.ITEM",
                EcrfFieldFunctions.name(studyEventOID, formOID, itemGroupOID, itemOID));
    }
}