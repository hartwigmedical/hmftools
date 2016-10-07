package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class ItemGroup implements OIDObject {

    @NotNull
    private final String OID;
    @NotNull
    private final List<String> itemOIDs;

    ItemGroup(@NotNull final String OID, @NotNull final List<String> itemOIDs) {
        this.OID = OID;
        this.itemOIDs = itemOIDs;
    }

    @NotNull
    public String OID() {
        return OID;
    }

    @NotNull
    List<String> itemOIDs() {
        return itemOIDs;
    }
}
