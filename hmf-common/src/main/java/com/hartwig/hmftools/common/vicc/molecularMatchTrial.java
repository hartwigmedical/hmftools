package com.hartwig.hmftools.common.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class molecularMatchTrial {

    public static StringBuilder readObjectMolecularMatchTrials(@NotNull JsonObject object) {
        //MolecularMatchTrials object
        StringBuilder stringToCSVMolecularMatchTrials = new StringBuilder();

        StringBuilder stringInterventionName = new StringBuilder();
        StringBuilder stringDescription = new StringBuilder();
        StringBuilder stringArmGroupLabel = new StringBuilder();
        StringBuilder stringInterventionType = new StringBuilder();
        StringBuilder stringOtherName = new StringBuilder();

        StringBuilder stringStatus = new StringBuilder();
        StringBuilder stringCity = new StringBuilder();
        StringBuilder stringValid = new StringBuilder();
        StringBuilder stringZip = new StringBuilder();
        StringBuilder stringCreated = new StringBuilder();
        StringBuilder stringCountry = new StringBuilder();
        StringBuilder stringId = new StringBuilder();
        StringBuilder stringLastUpdated = new StringBuilder();
        StringBuilder stringContact = new StringBuilder();
        StringBuilder stringState = new StringBuilder();
        StringBuilder stringStreet = new StringBuilder();
        StringBuilder stringLocation = new StringBuilder();
        StringBuilder stringPoBox = new StringBuilder();
        StringBuilder stringFailedGeocode = new StringBuilder();
        StringBuilder stringGeo = new StringBuilder();
        StringBuilder stringValidMessage = new StringBuilder();
        StringBuilder stringName = new StringBuilder();
        StringBuilder stringType = new StringBuilder();
        StringBuilder stringCoordinates = new StringBuilder();
        StringBuilder stringLat = new StringBuilder();
        StringBuilder stringLon = new StringBuilder();

        StringBuilder stringFacet = new StringBuilder();
        StringBuilder stringCompositeKey = new StringBuilder();
        StringBuilder stringSuppress = new StringBuilder();
        StringBuilder stringGenerateBy = new StringBuilder();
        StringBuilder stringFilterType = new StringBuilder();
        StringBuilder stringTerm = new StringBuilder();
        StringBuilder stringCustom = new StringBuilder();
        StringBuilder stringPriority = new StringBuilder();
        StringBuilder stringAlias = new StringBuilder();
        StringBuilder stringGeneratedByTerm = new StringBuilder();

        if (object.getAsJsonObject("molecularmatch_trials") != null) {
            List<String> keysOfMolecularMatchTrials = new ArrayList<>(object.getAsJsonObject("molecularmatch_trials").keySet());
            for (int x = 0; x < keysOfMolecularMatchTrials.size(); x++) {
                if (keysOfMolecularMatchTrials.get(x).equals("interventions")) {
                    JsonArray molecluarMatchTrialsArray =
                            object.getAsJsonObject("molecularmatch_trials").get(keysOfMolecularMatchTrials.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchTrialsArray.size(); i++) {
                        JsonObject objectinterventions = (JsonObject) molecluarMatchTrialsArray.get(i);
                        List<String> keysIntervations = new ArrayList<>(objectinterventions.keySet());
                        for (int p = 0; p < keysIntervations.size(); p++) {
                            if (keysIntervations.get(p).equals("intervention_name")) {
                                JsonElement interventionName = objectinterventions.get(keysIntervations.get(p));
                                stringInterventionName.append(interventionName).append(",");
                            } else if (keysIntervations.get(p).equals("other_name")) {
                                JsonElement otherName = objectinterventions.get(keysIntervations.get(p));
                                stringOtherName.append(otherName).append(",");
                            } else if (keysIntervations.get(p).equals("description")) {
                                JsonElement description = objectinterventions.get(keysIntervations.get(p));
                                stringDescription.append(description).append(",");
                            } else if (keysIntervations.get(p).equals("arm_group_label")) {
                                JsonElement armGroupLabel = objectinterventions.get(keysIntervations.get(p));
                                stringArmGroupLabel.append(armGroupLabel).append(",");
                            } else if (keysIntervations.get(p).equals("intervention_type")) {
                                JsonElement interventionType = objectinterventions.get(keysIntervations.get(p));
                                stringInterventionType.append(interventionType).append(",");
                            }
                        }
                    }
                    stringToCSVMolecularMatchTrials.append(stringInterventionName)
                            .append(";")
                            .append(stringOtherName)
                            .append(";")
                            .append(stringDescription)
                            .append(";")
                            .append(stringArmGroupLabel)
                            .append(";")
                            .append(stringInterventionType)
                            .append(";");
                } else if (keysOfMolecularMatchTrials.get(x).equals("locations")) {
                    JsonArray molecluarMatchTrialsArray =
                            object.getAsJsonObject("molecularmatch_trials").get(keysOfMolecularMatchTrials.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchTrialsArray.size(); i++) {
                        JsonObject objectLocations = (JsonObject) molecluarMatchTrialsArray.get(i);
                        List<String> keysLocations = new ArrayList<>(objectLocations.keySet());
                        for (int p = 0; p < keysLocations.size(); p++) {
                            if (keysLocations.get(p).equals("status")) {
                                JsonElement status = objectLocations.get(keysLocations.get(p));
                                stringStatus.append(status).append(",");
                            } else if (keysLocations.get(p).equals("city")) {
                                JsonElement city = objectLocations.get(keysLocations.get(p));
                                stringCity.append(city).append(",");
                            } else if (keysLocations.get(p).equals("_valid")) {
                                JsonElement valid = objectLocations.get(keysLocations.get(p));
                                stringValid.append(valid).append(",");
                            } else if (keysLocations.get(p).equals("zip")) {
                                JsonElement zip = objectLocations.get(keysLocations.get(p));
                                stringZip.append(zip).append(",");
                            } else if (keysLocations.get(p).equals("created")) {
                                JsonElement created = objectLocations.get(keysLocations.get(p));
                                stringCreated.append(created).append(",");
                            } else if (keysLocations.get(p).equals("id")) {
                                JsonElement id = objectLocations.get(keysLocations.get(p));
                                stringId.append(id).append(",");
                            } else if (keysLocations.get(p).equals("lastUpdated")) {
                                JsonElement lastUpdated = objectLocations.get(keysLocations.get(p));
                                stringLastUpdated.append(lastUpdated).append(",");
                            } else if (keysLocations.get(p).equals("contact")) {
                                JsonElement contact = objectLocations.get(keysLocations.get(p));
                                List<String> keysContact = Lists.newArrayList(contact.getAsJsonObject().keySet());
                                for (int r = 0; r < keysContact.size(); r++) {
                                    JsonElement contactData = contact.getAsJsonObject().get(keysContact.get(r));
                                    stringContact.append(contactData).append(";");
                                }
                            } else if (keysLocations.get(p).equals("state")) {
                                JsonElement state = objectLocations.get(keysLocations.get(p));
                                stringState.append(state).append(",");
                            } else if (keysLocations.get(p).equals("street")) {
                                JsonElement street = objectLocations.get(keysLocations.get(p));
                                stringStreet.append(street).append(",");
                            } else if (keysLocations.get(p).equals("location")) {
                                JsonElement location = objectLocations.get(keysLocations.get(p));
                                if (!location.isJsonNull()) {
                                    List<String> keys =
                                            new ArrayList<>(objectLocations.get(keysLocations.get(p)).getAsJsonObject().keySet());
                                    for (int y = 0; y < keys.size(); y++) {
                                        if (keys.get(y).equals("type")) {
                                            JsonElement type = location.getAsJsonObject().get(keys.get(y));
                                            stringType.append(type).append(",");
                                        } else if (keys.get(y).equals("coordinates")) {
                                            JsonElement coordinates = location.getAsJsonObject().get(keys.get(y));
                                            stringCoordinates.append(coordinates).append(",");
                                        }
                                    }
                                }
                            } else if (keysLocations.get(p).equals("po_box")) {
                                JsonElement po_box = objectLocations.get(keysLocations.get(p));
                                stringPoBox.append(po_box).append(",");
                            } else if (keysLocations.get(p).equals("failedGeocode")) {
                                JsonElement failedGeocode = objectLocations.get(keysLocations.get(p));
                                stringFailedGeocode.append(failedGeocode).append(",");
                            } else if (keysLocations.get(p).equals("geo")) {
                                JsonElement geo = objectLocations.get(keysLocations.get(p));
                                List<String> keys = new ArrayList<>(objectLocations.get(keysLocations.get(p)).getAsJsonObject().keySet());
                                for (int v = 0; v < keys.size(); v++) {
                                    if (keys.get(v).equals("lat")) {
                                        JsonElement lat = geo.getAsJsonObject().get(keys.get(v));
                                        stringLat.append(lat).append(",");
                                    } else if (keys.get(v).equals("lon")) {
                                        JsonElement lon = geo.getAsJsonObject().get(keys.get(v));
                                        stringLon.append(lon).append(",");
                                    }
                                }
                                stringGeo.append(geo).append(",");
                            } else if (keysLocations.get(p).equals("_validMessage")) {
                                JsonElement validMessage = objectLocations.get(keysLocations.get(p));
                                stringValidMessage.append(validMessage).append(",");
                            } else if (keysLocations.get(p).equals("name")) {
                                JsonElement name = objectLocations.get(keysLocations.get(p));
                                stringName.append(name).append(",");
                            }
                        }
                    }
                    stringToCSVMolecularMatchTrials.append(stringStatus)
                            .append(";")
                            .append(stringCity)
                            .append(";")
                            .append(stringValid)
                            .append(";")
                            .append(stringZip)
                            .append(";")
                            .append(stringCreated)
                            .append(";")
                            .append(stringCountry)
                            .append(";")
                            .append(stringId)
                            .append(";")
                            .append(stringLastUpdated)
                            .append(";")
                            .append(stringContact)
                            .append(stringState)
                            .append(";")
                            .append(stringStreet)
                            .append(";")
                            .append(stringType)
                            .append(";")
                            .append(stringCoordinates)
                            .append(";")
                            .append(stringPoBox)
                            .append(";")
                            .append(stringFailedGeocode)
                            .append(";")
                            .append(stringLat)
                            .append(";")
                            .append(stringLon)
                            .append(";")
                            .append(stringValidMessage)
                            .append(";")
                            .append(stringName)
                            .append(";");
                } else if (keysOfMolecularMatchTrials.get(x).equals("overallContact")) {
                    if (!object.getAsJsonObject("molecularmatch_trials").get(keysOfMolecularMatchTrials.get(x)).isJsonNull()) {
                        JsonObject objectOverallContact =
                                object.getAsJsonObject("molecularmatch_trials").get(keysOfMolecularMatchTrials.get(x)).getAsJsonObject();
                        List<String> keysOfOverallContact = new ArrayList<>(objectOverallContact.keySet());
                        for (int u = 0; u < keysOfOverallContact.size(); u++) {
                            if (keysOfOverallContact.get(u).equals("phone")) {
                                stringToCSVMolecularMatchTrials.append(objectOverallContact.get(keysOfOverallContact.get(u))).append(";");
                            } else if (keysOfOverallContact.get(u).equals("last_name")) {
                                if (u == 0) {
                                    stringToCSVMolecularMatchTrials.append(";");
                                }
                                stringToCSVMolecularMatchTrials.append(objectOverallContact.get(keysOfOverallContact.get(u))).append(";");
                            } else if (keysOfOverallContact.get(u).equals("email")) {
                                stringToCSVMolecularMatchTrials.append(objectOverallContact.get(keysOfOverallContact.get(u))).append(";");
                            } else if (keysOfOverallContact.get(u).equals("affiliation")) {
                                stringToCSVMolecularMatchTrials.append(objectOverallContact.get(keysOfOverallContact.get(u))).append(";");
                            } else {
                                stringToCSVMolecularMatchTrials.append(";;;;");
                            }
                        }
                    }
                } else if (keysOfMolecularMatchTrials.get(x).equals("tags")) {
                    JsonArray molecluarMatchTrialsArray =
                            object.getAsJsonObject("molecularmatch_trials").get(keysOfMolecularMatchTrials.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchTrialsArray.size(); i++) {
                        JsonObject objectTags = (JsonObject) molecluarMatchTrialsArray.get(i);
                        List<String> keysTags = new ArrayList<>(objectTags.keySet());
                        for (int p = 0; p < keysTags.size(); p++) {
                            if (keysTags.get(p).equals("facet")) {
                                JsonElement facet = objectTags.get(keysTags.get(p));
                                stringFacet.append(facet).append(",");
                            } else if (keysTags.get(p).equals("compositeKey")) {
                                JsonElement compositeKey = objectTags.get(keysTags.get(p));
                                stringCompositeKey.append(compositeKey).append(",");
                            } else if (keysTags.get(p).equals("suppress")) {
                                JsonElement suppress = objectTags.get(keysTags.get(p));
                                stringSuppress.append(suppress).append(",");
                            } else if (keysTags.get(p).equals("generateBy")) {
                                JsonElement generateBy = objectTags.get(keysTags.get(p));
                                stringGenerateBy.append(generateBy).append(",");
                            } else if (keysTags.get(p).equals("filterType")) {
                                JsonElement filterType = objectTags.get(keysTags.get(p));
                                stringFilterType.append(filterType).append(",");
                            } else if (keysTags.get(p).equals("term")) {
                                JsonElement term = objectTags.get(keysTags.get(p));
                                stringTerm.append(term).append(",");
                            } else if (keysTags.get(p).equals("custom")) {
                                JsonElement custom = objectTags.get(keysTags.get(p));
                                stringCustom.append(custom).append(",");
                            } else if (keysTags.get(p).equals("priority")) {
                                JsonElement priority = objectTags.get(keysTags.get(p));
                                stringPriority.append(priority).append(",");
                            } else if (keysTags.get(p).equals("alias")) {
                                JsonElement alias = objectTags.get(keysTags.get(p));
                                stringAlias.append(alias).append(",");
                            } else if (keysTags.get(p).equals("generatedByTerm")) {
                                JsonElement generatedByTerm = objectTags.get(keysTags.get(p));
                                stringGeneratedByTerm.append(generatedByTerm).append(",");
                            }
                        }
                    }
                    stringToCSVMolecularMatchTrials
                            .append(stringFacet)
                            .append(";")
                            .append(stringCompositeKey)
                            .append(";")
                            .append(stringSuppress)
                            .append(";")
                            .append(stringFilterType)
                            .append(";")
                            .append(stringTerm)
                            .append(";")
                            .append(stringCustom)
                            .append(";")
                            .append(stringPriority)
                            .append(";")
                            .append(stringAlias)
                            .append(";")
                            .append(stringGeneratedByTerm)
                            .append(";")
                            .append(stringGenerateBy)
                            .append(";");
                } else {
                    stringToCSVMolecularMatchTrials.append(object.getAsJsonObject("molecularmatch_trials")
                            .get(keysOfMolecularMatchTrials.get(x))).append(";");
                }
            }
        }
        return stringToCSVMolecularMatchTrials;
    }
}
