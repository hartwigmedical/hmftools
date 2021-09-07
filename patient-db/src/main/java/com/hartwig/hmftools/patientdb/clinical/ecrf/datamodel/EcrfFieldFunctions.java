package com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel;

import org.jetbrains.annotations.NotNull;

final class EcrfFieldFunctions {

    private static final String OID_SEPARATOR = ".";

    private EcrfFieldFunctions() {
    }

    @NotNull
    static String name(@NotNull String studyEventOID, @NotNull String formOID, @NotNull String itemGroupOID,
            @NotNull String itemOID) {
        String study = lastElement(studyEventOID);
        String form = lastElement(formOID);
        String itemGroup = lastElement(itemGroupOID);
        String item = lastElement(itemOID);

        return (study + OID_SEPARATOR + form + OID_SEPARATOR + itemGroup + OID_SEPARATOR + item).toUpperCase();
    }

    @NotNull
    private static String lastElement(@NotNull String OID) {
        String[] elements = OID.split("\\" + OID_SEPARATOR);
        return elements[elements.length - 1];
    }
}
