package com.hartwig.hmftools.common.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class brca {

    public static StringBuilder readObjectBRCA(@NotNull JsonObject object) {
        //brca object
        StringBuilder stringToCSVBRCA = new StringBuilder();
        if (object.getAsJsonObject("brca") != null) {
            for (int i = 0; i < object.getAsJsonObject("brca").keySet().size(); i++) {
                List<String> keysOfBRCAObject = new ArrayList<>(object.getAsJsonObject("brca").keySet());
                stringToCSVBRCA.append(object.getAsJsonObject("brca").get(keysOfBRCAObject.get(i))).append(";"); // brca data
            }
        } else {
            stringToCSVBRCA.append(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;"
                    + ";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;");
        }

        return stringToCSVBRCA;
    }
}
