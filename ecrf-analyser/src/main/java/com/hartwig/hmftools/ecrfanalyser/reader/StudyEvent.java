package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class StudyEvent {

    @NotNull
    private final String OID;
    @NotNull
    private final List<String> formOIDs;

    StudyEvent(@NotNull final String OID, @NotNull final List<String> formOIDs) {
        this.OID = OID;
        this.formOIDs = formOIDs;
    }
}
