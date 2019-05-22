package com.hartwig.hmftools.common.vicc;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class tags {

    public static StringBuilder readObjectTags(@NotNull JsonObject object) {
        //Tags object
        StringBuilder stringToCSVTags = new StringBuilder();
        JsonArray arrayTags = object.getAsJsonArray("tags");
        String tags = arrayTags.toString();
        tags = tags.substring(1, tags.length() - 1);
        stringToCSVTags.append(tags).append(";"); // tags data
        return stringToCSVTags;
    }
}
