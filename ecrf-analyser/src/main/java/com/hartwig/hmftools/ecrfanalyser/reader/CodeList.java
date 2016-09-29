package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class CodeList {

    @NotNull
    private final String OID;
    @NotNull
    private final List<String> codeListItems;

    public CodeList(@NotNull final String OID, @NotNull final List<String> codeListItems) {
        this.OID = OID;
        this.codeListItems = codeListItems;
    }

    @Override
    public String toString() {
        return "CodeList{" + "OID='" + OID + '\'' + ", codeListItems=" + codeListItems + '}';
    }
}
