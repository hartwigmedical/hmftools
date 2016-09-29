package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.Map;

import org.jetbrains.annotations.NotNull;

class CodeList {

    @NotNull
    private final String OID;
    @NotNull
    private final Map<Integer, String> values;

    CodeList(@NotNull final String OID, @NotNull final Map<Integer, String> values) {
        this.OID = OID;
        this.values = values;
    }

    @NotNull
    String OID() {
        return OID;
    }

    @NotNull
    Map<Integer, String> values() {
        return values;
    }

    @Override
    public String toString() {
        return "CodeList{" + "OID='" + OID + '\'' + ", values=" + values + '}';
    }
}
