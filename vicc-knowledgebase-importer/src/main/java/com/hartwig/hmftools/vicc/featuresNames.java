package com.hartwig.hmftools.vicc;

import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class featuresNames {

    public static void readObjectFeaturesNames(@NotNull JsonObject object) {
        //feature_names object

        if (object.get("feature_names") == null) {
            //TODO: add empty field to sql
        } else {
            JsonElement element = object.get("feature_names");
            //TODO: add to sql
        }
    }
}
