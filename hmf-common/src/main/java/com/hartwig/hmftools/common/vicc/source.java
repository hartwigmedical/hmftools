package com.hartwig.hmftools.common.vicc;

import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;

import org.jetbrains.annotations.NotNull;

public class source {

    public static void readObjectSource(@NotNull JsonObject object) {
        //Source object
        JsonPrimitive source = object.getAsJsonPrimitive("source"); // source data
        // TODO: add sql
    }
}
