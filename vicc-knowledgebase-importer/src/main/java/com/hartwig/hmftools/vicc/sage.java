package com.hartwig.hmftools.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class sage {


    public static void readObjectSage(@NotNull JsonObject object) {
        //SAGE object
        List<String> keysOfSage;

        if (object.getAsJsonObject("sage") != null) {
            keysOfSage = new ArrayList<>(object.getAsJsonObject("sage").keySet());

            for (int j = 0; j < keysOfSage.size(); j++) {
                object.getAsJsonObject("sage").get(keysOfSage.get(j));
                // TODO: add to sql
            }
        }
    }
}
