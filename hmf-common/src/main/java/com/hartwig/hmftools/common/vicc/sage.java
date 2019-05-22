package com.hartwig.hmftools.common.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class sage {

    public static StringBuilder readObjectSage(@NotNull JsonObject object) {
        //SAGE object
        StringBuilder stringToCSVSage = new StringBuilder();
        List<String> keysOfSage;

        if (object.getAsJsonObject("sage") != null) {
            keysOfSage = new ArrayList<>(object.getAsJsonObject("sage").keySet());

            for (int j = 0; j < keysOfSage.size(); j++) {
                stringToCSVSage.append(object.getAsJsonObject("sage").get(keysOfSage.get(j))).append(";");
            }
        } else {
            stringToCSVSage.append(";;;;;;;;");
        }
        return stringToCSVSage;
    }
}
