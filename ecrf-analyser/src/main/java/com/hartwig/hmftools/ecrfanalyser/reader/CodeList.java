package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;
import java.util.Map;

import org.jetbrains.annotations.NotNull;

class CodeList {

    @NotNull
    private final String OID;
    @NotNull
    private final Map<Integer, String> codeListItems;

    CodeList(@NotNull final String OID, @NotNull final Map<Integer, String> codeListItems) {
        this.OID = OID;
        this.codeListItems = codeListItems;
    }

    @Override
    public String toString() {
        return "CodeList{" + "OID='" + OID + '\'' + ", codeListItems=" + codeListItems + '}';
    }
}
