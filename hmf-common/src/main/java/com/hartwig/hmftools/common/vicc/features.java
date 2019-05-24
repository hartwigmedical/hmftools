package com.hartwig.hmftools.common.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class features {

    public static StringBuilder readObjectFeatures(@NotNull JsonObject object) {
        //features
        StringBuilder stringToCSVFeatures = new StringBuilder();
        JsonArray arrayFeatures = object.get("features").getAsJsonArray();
        for (int p = 0; p < arrayFeatures.size(); p++) {
            JsonObject objectFeatures = (JsonObject) arrayFeatures.get(p);
            List<String> keysFeatures = new ArrayList<>(objectFeatures.keySet());
            for (int i = 0; i < keysFeatures.size(); i++) {
                if (keysFeatures.get(i).equals("sequence_ontology")) {
                    JsonObject objectSequenceOntolgy = objectFeatures.get(keysFeatures.get(i)).getAsJsonObject();
                    List<String> keysSequenceOntology = new ArrayList<>(objectSequenceOntolgy.keySet());
                    for (int z = 0; z < keysSequenceOntology.size(); z++) {
                        stringToCSVFeatures.append(objectSequenceOntolgy.get(keysSequenceOntology.get(z))).append(";");
                    }

                } else {
                    stringToCSVFeatures.append(objectFeatures.get(keysFeatures.get(i))).append(";");
                }

            }
        }
        return stringToCSVFeatures;
    }
}
