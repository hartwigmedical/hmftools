package com.hartwig.hmftools.common.vicc;

import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class source {

    public static StringBuilder readObjectSource(@NotNull JsonObject object) {
        //Source object
        StringBuilder stringToCSVSource = new StringBuilder();
        stringToCSVSource.append(object.getAsJsonPrimitive("source")).append(";"); // source data
        return stringToCSVSource;
    }
}
