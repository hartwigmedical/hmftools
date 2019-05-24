package com.hartwig.hmftools.common.vicc;

import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class featuresNames {

    public static StringBuilder readObjectFeaturesNames(@NotNull JsonObject object) {
        //feature_names object
        StringBuilder stringToCSVFeaturesNames = new StringBuilder();

        if (object.get("feature_names") == null) {
            stringToCSVFeaturesNames.append("null").append(";");
        } else {
            String element = object.get("feature_names").toString().replaceAll(";", ":");
            stringToCSVFeaturesNames.append(element).append(";");
        }
        return stringToCSVFeaturesNames;
    }
}
