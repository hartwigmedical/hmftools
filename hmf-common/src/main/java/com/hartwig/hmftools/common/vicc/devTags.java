package com.hartwig.hmftools.common.vicc;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class devTags {

    public static StringBuilder readObjectDevTags(@NotNull JsonObject object) {
        //dev_tags
        StringBuilder stringToCSVDevTags = new StringBuilder();
        JsonArray arrayDevTags = object.getAsJsonArray("dev_tags");
        String devTags = arrayDevTags.toString();
        devTags = devTags.substring(1, devTags.length() - 1);
        stringToCSVDevTags.append(devTags).append(";"); // dev tags data
        return stringToCSVDevTags;
    }
}
