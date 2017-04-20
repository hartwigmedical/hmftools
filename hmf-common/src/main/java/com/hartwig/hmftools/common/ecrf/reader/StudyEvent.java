package com.hartwig.hmftools.common.ecrf.reader;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class StudyEvent implements OIDObject {

    @NotNull
    private final String OID;
    @NotNull
    private final List<String> formOIDs;

    StudyEvent(@NotNull final String OID, @NotNull final List<String> formOIDs) {
        this.OID = OID;
        this.formOIDs = formOIDs;
    }

    @NotNull
    public String OID() {
        return OID;
    }

    @NotNull
    List<String> formOIDs() {
        return formOIDs;
    }
}
