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

    private static final Logger LOGGER = LogManager.getLogger(sage.class);

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

    public static StringBuilder readObjectSageSpecificFields(@NotNull JsonObject object) {
        //SAGE object
        StringBuilder stringToCSVSage = new StringBuilder();
        List<String> keysOfSage;
        List<String> listSage = Lists.newArrayList();

        if (object.getAsJsonObject("sage") != null) {
            keysOfSage = new ArrayList<>(object.getAsJsonObject("sage").keySet());

            for (int j = 0; j < keysOfSage.size(); j++) {
                if (keysOfSage.get(j).equals("publication_url")) {
                    listSage.add(0, object.getAsJsonObject("sage").get(keysOfSage.get(j)).toString());
                }
                if (keysOfSage.get(j).equals("clinical_manifestation")) {
                    listSage.add(0, object.getAsJsonObject("sage").get(keysOfSage.get(j)).toString());
                }
                if (keysOfSage.get(j).equals("drug_labels")) {
                    listSage.add(0, object.getAsJsonObject("sage").get(keysOfSage.get(j)).toString());
                }
                if (keysOfSage.get(j).equals("evidence_label")) {
                    listSage.add(0, object.getAsJsonObject("sage").get(keysOfSage.get(j)).toString());
                }
                if (keysOfSage.get(j).equals("response_type")) {
                    listSage.add(0, object.getAsJsonObject("sage").get(keysOfSage.get(j)).toString());
                }
                if (keysOfSage.get(j).equals("gene")) {
                    listSage.add(0, object.getAsJsonObject("sage").get(keysOfSage.get(j)).toString());
                }
            }
               stringToCSVSage.append(String.join(";" , listSage));
        }
        LOGGER.info(stringToCSVSage);
        return stringToCSVSage;
    }
}
