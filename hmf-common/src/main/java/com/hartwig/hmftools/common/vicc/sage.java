package com.hartwig.hmftools.common.vicc;

import java.security.acl.LastOwnerException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
