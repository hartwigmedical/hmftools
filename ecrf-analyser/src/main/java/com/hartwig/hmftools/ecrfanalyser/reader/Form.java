package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class Form implements OIDObject {

    @NotNull
    private final String OID;
    @NotNull
    private final List<String> itemGroupOIDs;

    Form(@NotNull final String OID, @NotNull final List<String> itemGroupOIDs) {
        this.OID = OID;
        this.itemGroupOIDs = itemGroupOIDs;
    }

    @NotNull
    public String OID() {
        return OID;
    }

    @NotNull
    List<String> itemGroupOIDs() {
        return itemGroupOIDs;
    }
}
