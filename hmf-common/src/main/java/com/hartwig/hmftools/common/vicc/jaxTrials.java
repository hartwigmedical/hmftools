package com.hartwig.hmftools.common.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class jaxTrials {

    public static StringBuilder readObjectJaxTrials(@NotNull JsonObject object) {
        //Jax_trails object
        StringBuilder stringToCSVJaxTrials = new StringBuilder();
        StringBuilder stringSource = new StringBuilder();
        StringBuilder stringId = new StringBuilder();
        StringBuilder stringName = new StringBuilder();
        StringBuilder stringTherapyName = new StringBuilder();
        StringBuilder stringTherapyId = new StringBuilder();
        StringBuilder stringRequirementType = new StringBuilder();
        StringBuilder stringProfileId = new StringBuilder();
        StringBuilder stringProfileName = new StringBuilder();

        if (object.getAsJsonObject("jax_trials") != null) {
            List<String> keysOfJaxTrials = new ArrayList<>(object.getAsJsonObject("jax_trials").keySet());
            for (int x = 0; x < keysOfJaxTrials.size(); x++) {
                if (keysOfJaxTrials.get(x).equals("indications")) {
                    JsonArray jaxTrailsArray = object.getAsJsonObject("jax_trials").get(keysOfJaxTrials.get(x)).getAsJsonArray();
                    for (int u = 0; u < jaxTrailsArray.size(); u++) {
                        JsonObject objectIndications = (JsonObject) jaxTrailsArray.get(u);
                        List<String> keysIndications = new ArrayList<>(objectIndications.keySet());
                        for (int v = 0; v < objectIndications.keySet().size(); v++) {
                            if (keysIndications.get(v).equals("source")) {
                                JsonElement source = objectIndications.get(keysIndications.get(v));
                                stringSource.append(source).append(",");
                            } else if (keysIndications.get(v).equals("id")) {
                                JsonElement id = objectIndications.get(keysIndications.get(v));
                                stringId.append(id).append(",");
                            } else if (keysIndications.get(v).equals("name")) {
                                JsonElement name = objectIndications.get(keysIndications.get(v));
                                stringName.append(name).append(",");
                            }
                        }
                    }
                    stringToCSVJaxTrials.append(stringSource).append(";").append(stringId).append(";").append(stringName).append(";");
                } else if (keysOfJaxTrials.get(x).equals("variantRequirementDetails")) {
                    JsonArray jaxTrailsArray = object.getAsJsonObject("jax_trials").get(keysOfJaxTrials.get(x)).getAsJsonArray();
                    for (int u = 0; u < jaxTrailsArray.size(); u++) {
                        JsonObject objectVariantRequirementDetails = (JsonObject) jaxTrailsArray.get(u);
                        List<String> keysVariantRequirementDetails = new ArrayList<>(objectVariantRequirementDetails.keySet());
                        for (int v = 0; v < keysVariantRequirementDetails.size(); v++) {
                            if (keysVariantRequirementDetails.get(v).equals("molecularProfile")) {
                                JsonObject objectMolecularProfile =
                                        objectVariantRequirementDetails.get(keysVariantRequirementDetails.get(v)).getAsJsonObject();
                                List<String> keysMolecularProfile = new ArrayList<>(objectMolecularProfile.keySet());
                                for (int p = 0; p < keysMolecularProfile.size(); p++) {
                                    if (keysMolecularProfile.get(p).equals("profileName")) {
                                        JsonElement profileName = objectMolecularProfile.get(keysMolecularProfile.get(p));
                                        stringProfileName.append(profileName).append(",");
                                    } else if (keysMolecularProfile.get(p).equals("id")) {
                                        JsonElement id = objectMolecularProfile.get(keysMolecularProfile.get(p));
                                        stringId.append(id).append(",");
                                    }
                                }
                            } else if (keysVariantRequirementDetails.get(v).equals("requirementType")) {
                                JsonElement requireType = objectVariantRequirementDetails.get(keysVariantRequirementDetails.get(v));
                                stringRequirementType.append(requireType).append(",");
                            }
                        }
                    }
                    stringToCSVJaxTrials.append(stringProfileName)
                            .append(";")
                            .append(stringProfileId)
                            .append(";")
                            .append(stringRequirementType)
                            .append(";");
                } else if (keysOfJaxTrials.get(x).equals("therapies")) {
                    JsonElement jaxTrailsArray = object.getAsJsonObject("jax_trials").get(keysOfJaxTrials.get(x));
                    for (int v = 0; v < jaxTrailsArray.getAsJsonArray().size(); v++) {
                        List<String> keysTherapies = new ArrayList<>(jaxTrailsArray.getAsJsonArray().get(v).getAsJsonObject().keySet());
                        JsonObject objectTherapies = jaxTrailsArray.getAsJsonArray().get(v).getAsJsonObject();
                        for (int z = 0; z < keysTherapies.size(); z++) {
                            if (keysTherapies.get(z).equals("id")) {
                                JsonElement id = objectTherapies.get(keysTherapies.get(z));
                                stringTherapyId.append(id).append(",");
                            } else if (keysTherapies.get(z).equals("therapyName")) {
                                JsonElement therapyName = objectTherapies.get(keysTherapies.get(z));
                                stringTherapyName.append(therapyName).append(",");
                            }
                        }
                    }
                    stringToCSVJaxTrials.append(stringTherapyId).append(";").append(stringTherapyName).append(";");
                } else {
                    stringToCSVJaxTrials.append(object.getAsJsonObject("jax_trials").get(keysOfJaxTrials.get(x))).append(";");
                }
            }
        } else {
            stringToCSVJaxTrials.append(";;;;;;;;;;;;;;;;");
        }
        return stringToCSVJaxTrials;
    }
}
