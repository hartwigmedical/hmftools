package com.hartwig.hmftools.common.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class pmkb {

    public static StringBuilder readObjectPmkb(@NotNull JsonObject object) {
        //PMKB object
        StringBuilder stringToCSVPmkb = new StringBuilder();
        List<String> keysOfPmkb;

        StringBuilder stringId = new StringBuilder();
        StringBuilder stringName = new StringBuilder();
        if (object.getAsJsonObject("pmkb") != null) {
            keysOfPmkb = new ArrayList<>(object.getAsJsonObject("pmkb").keySet());
            for (int j = 0; j < keysOfPmkb.size(); j++) {
                if (keysOfPmkb.get(j).equals("tumor")) {
                    JsonObject pmkbObject = object.getAsJsonObject("pmkb").get(keysOfPmkb.get(j)).getAsJsonObject();
                    List<String> keysOfVariant = new ArrayList<>(pmkbObject.keySet());
                    for (int x = 0; x < pmkbObject.keySet().size(); x++) {
                        stringToCSVPmkb.append(pmkbObject.get(keysOfVariant.get(x))).append(";");
                    }
                } else if (keysOfPmkb.get(j).equals("tissues")) {
                    JsonArray arrayTissue = object.getAsJsonObject("pmkb").get(keysOfPmkb.get(j)).getAsJsonArray();
                    for (int x = 0; x < arrayTissue.size(); x++) {
                        JsonObject objectTissue = (JsonObject) arrayTissue.get(x);

                        for (int v = 0; v < objectTissue.keySet().size(); v++) {
                            List<String> keysTissue = new ArrayList<>(objectTissue.keySet());
                            if (keysTissue.get(v).equals("id")) {
                                JsonElement idTissue = objectTissue.get(keysTissue.get(v));
                                stringId.append(idTissue).append(",");
                            } else if (keysTissue.get(v).equals("name")) {
                                JsonElement nameTissue = objectTissue.get(keysTissue.get(v));
                                stringName.append(nameTissue).append(",");
                            }
                        }
                    }
                    stringToCSVPmkb.append(stringId).append(";").append(stringName).append(";");
                } else if (keysOfPmkb.get(j).equals("variant")) {
                    JsonObject pmkbObject = object.getAsJsonObject("pmkb").get(keysOfPmkb.get(j)).getAsJsonObject();
                    List<String> keysOfVariant = new ArrayList<>(pmkbObject.keySet());
                    for (int x = 0; x < pmkbObject.keySet().size(); x++) {
                        if (keysOfVariant.get(x).equals("gene")) {
                            JsonElement elementGene = object.getAsJsonObject("pmkb").get("variant");
                            List<String> keysGene = new ArrayList<>(elementGene.getAsJsonObject().get("gene").getAsJsonObject().keySet());

                            for (int d = 0; d < keysGene.size(); d++) {
                                stringToCSVPmkb.append(elementGene.getAsJsonObject() // association data
                                        .get("gene").getAsJsonObject().get(keysGene.get(d))).append(";");
                            }
                        } else {
                            stringToCSVPmkb.append(pmkbObject.get(keysOfVariant.get(x))).append(";");
                        }
                    }
                }
            }
        } else {
            stringToCSVPmkb.append(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;");
        }

        return stringToCSVPmkb;
    }
}
