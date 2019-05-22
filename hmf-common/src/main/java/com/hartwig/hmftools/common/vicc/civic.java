package com.hartwig.hmftools.common.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class civic {

    public static StringBuilder readObjectCIVIC(@NotNull JsonObject object) {
        //CIVIC object
        StringBuilder stringToCSVCIVIC = new StringBuilder();
        StringBuilder stringDisplayName = new StringBuilder();
        StringBuilder stringDescription = new StringBuilder();
        StringBuilder stringUrl = new StringBuilder();
        StringBuilder stringSoId = new StringBuilder();
        StringBuilder stringId = new StringBuilder();
        StringBuilder stringName = new StringBuilder();

        if (object.getAsJsonObject("civic") != null) {
            List<String> keysOfCivic = new ArrayList<>(object.getAsJsonObject("civic").keySet());
            for (int x = 0; x < keysOfCivic.size(); x++) {
                if (keysOfCivic.get(x).equals("variant_types")) {
                    JsonArray civicVariantTypesArray = object.getAsJsonObject("civic").get(keysOfCivic.get(x)).getAsJsonArray();
                    for (int z = 0; z < civicVariantTypesArray.size(); z++) {
                        JsonObject objectVariantTypes = (JsonObject) civicVariantTypesArray.get(z);
                        List<String> keysVariantTypes = new ArrayList<>(objectVariantTypes.keySet());
                        for (int i = 0; i < keysVariantTypes.size(); i++) {
                            if (keysVariantTypes.get(i).equals("display_name")) {
                                JsonElement displayName = objectVariantTypes.get(keysVariantTypes.get(i));
                                stringDisplayName.append(displayName).append(",");
                            } else if (keysVariantTypes.get(i).equals("description")) {
                                JsonElement discription = objectVariantTypes.get(keysVariantTypes.get(i));
                                stringDescription.append(discription).append(",");
                            } else if (keysVariantTypes.get(i).equals("url")) {
                                JsonElement url = objectVariantTypes.get(keysVariantTypes.get(i));
                                stringUrl.append(url).append(",");
                            } else if (keysVariantTypes.get(i).equals("so_id")) {
                                JsonElement soId = objectVariantTypes.get(keysVariantTypes.get(i));
                                stringSoId.append(soId).append(",");
                            } else if (keysVariantTypes.get(i).equals("id")) {
                                JsonElement id = objectVariantTypes.get(keysVariantTypes.get(i));
                                stringId.append(id).append(",");
                            } else if (keysVariantTypes.get(i).equals("name")) {
                                JsonElement name = objectVariantTypes.get(keysVariantTypes.get(i));
                                stringName.append(name).append(",");
                            }
                        }
                    }
                    stringToCSVCIVIC.append(stringDisplayName)
                            .append(";")
                            .append(stringDescription)
                            .append(";")
                            .append(stringUrl)
                            .append(";")
                            .append(stringSoId)
                            .append(";")
                            .append(stringId)
                            .append(";")
                            .append(stringName)
                            .append(";");
                } else if (keysOfCivic.get(x).equals("lifecycle_actions")) {
                    JsonObject civicObject = object.getAsJsonObject("civic").get(keysOfCivic.get(x)).getAsJsonObject();
                    List<String> keysOfLifeCycleActions = new ArrayList<>(civicObject.keySet());
                    for (int i = 0; i < keysOfLifeCycleActions.size(); i++) {
                        if (keysOfLifeCycleActions.get(i).equals("last_commented_on")) {
                            commandCIVIC(civicObject, "last_commented_on", stringToCSVCIVIC);
                        }
                        if (keysOfLifeCycleActions.get(i).equals("last_modified")) {
                            commandCIVIC(civicObject, "last_modified", stringToCSVCIVIC);
                        }
                        if (keysOfLifeCycleActions.get(i).equals("last_reviewed")) {
                            commandCIVIC(civicObject, "last_reviewed", stringToCSVCIVIC);
                        }
                    }

                } else if (keysOfCivic.get(x).equals("evidence_items")) {
                    JsonArray civicEvidenceItemsArray = object.getAsJsonObject("civic").get(keysOfCivic.get(x)).getAsJsonArray();
                    for (int z = 0; z < civicEvidenceItemsArray.size(); z++) {
                        JsonObject objectEvidenceitems = (JsonObject) civicEvidenceItemsArray.get(z);
                        List<String> keysEvidenceItems = new ArrayList<>(objectEvidenceitems.keySet());
                        for (int h = 0; h < keysEvidenceItems.size(); h++) {
                            if (keysEvidenceItems.get(h).equals("drugs")) {
                                JsonArray ArrayDrugs = objectEvidenceitems.get(keysEvidenceItems.get(h)).getAsJsonArray();
                                for (int r = 0; r < ArrayDrugs.size(); r++) {
                                    JsonObject objectDrugs = (JsonObject) ArrayDrugs.get(r);
                                    List<String> keysDrugs = new ArrayList<>(objectDrugs.keySet());
                                    for (int f = 0; f < keysDrugs.size(); f++) {
                                        stringToCSVCIVIC.append(objectDrugs.get(keysDrugs.get(f))).append(";");
                                    }
                                }
                            } else if (keysEvidenceItems.get(h).equals("disease")) {
                                JsonObject objectDiseases = objectEvidenceitems.get(keysEvidenceItems.get(h)).getAsJsonObject();
                                List<String> keysDiseases = new ArrayList<>(objectDiseases.keySet());
                                for (int a = 0; a < keysDiseases.size(); a++) {
                                    stringToCSVCIVIC.append(objectDiseases.get(keysDiseases.get(a))).append(";");
                                }
                            } else if (keysEvidenceItems.get(h).equals("source")) {
                                JsonObject objectSource = objectEvidenceitems.get(keysEvidenceItems.get(h)).getAsJsonObject();
                                List<String> keysSource = new ArrayList<>(objectSource.keySet());
                                for (int o = 0; o < keysSource.size(); o++) {
                                    if (keysSource.get(o).equals("publication_date")) {
                                        JsonObject objectPublication = objectSource.get(keysSource.get(o)).getAsJsonObject();
                                        List<String> keysPublication = new ArrayList<>(objectPublication.keySet());
                                        for (int d = 0; d < keysPublication.size(); d++) {
                                            stringToCSVCIVIC.append(objectPublication.get(keysPublication.get(d))).append(";");
                                        }
                                    } else {
                                        stringToCSVCIVIC.append(objectSource.get(keysSource.get(o))).append(";");
                                    }
                                }
                            } else {
                                stringToCSVCIVIC.append(objectEvidenceitems.get(keysEvidenceItems.get(h))).append(";");
                            }
                        }
                    }
                } else if (keysOfCivic.get(x).equals("coordinates")) {
                    JsonObject civicObject = object.getAsJsonObject("civic").get(keysOfCivic.get(x)).getAsJsonObject();
                    List<String> keysCoordinates = new ArrayList<>(civicObject.keySet());
                    for (int w = 0; w < keysCoordinates.size(); w++) {
                        stringToCSVCIVIC.append(civicObject.get(keysCoordinates.get(w))).append(";");

                    }
                } else {
                    stringToCSVCIVIC.append(object.getAsJsonObject("civic").get(keysOfCivic.get(x))).append(";");
                }
            }
        } else {
            stringToCSVCIVIC.append(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; "
                    + ";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;");
        }
        return stringToCSVCIVIC;
    }

    private static void commandCIVIC(@NotNull JsonObject civicObject, @NotNull String i, @NotNull StringBuilder stringToCSVCIVIC) {
        JsonObject objectLastCommened = civicObject.get(i).getAsJsonObject();
        List<String> keysLastCommened = new ArrayList<>(objectLastCommened.keySet());
        for (int p = 0; p < keysLastCommened.size(); p++) {
            if (keysLastCommened.get(p).equals("timestamp")) {
                stringToCSVCIVIC.append(objectLastCommened.get(i)).append(";");
            } else if (keysLastCommened.get(p).equals("user")) {
                JsonObject objectUser = objectLastCommened.get(keysLastCommened.get(p)).getAsJsonObject();
                List<String> keysUser = new ArrayList<>(objectUser.keySet());
                for (int g = 0; g < keysUser.size(); g++) {
                    if (keysUser.get(g).equals("organization")) {
                        JsonObject objectOrganization = objectUser.get(keysUser.get(g)).getAsJsonObject();
                        List<String> keysOrganization = new ArrayList<>(objectOrganization.keySet());
                        for (int d = 0; d < keysOrganization.size(); d++) {
                            if (keysOrganization.get(d).equals("profile_image")) {
                                JsonObject objectProfile_image = objectOrganization.get(keysOrganization.get(d)).getAsJsonObject();
                                List<String> keysProfileImage = new ArrayList<>(objectProfile_image.keySet());
                                for (int t = 0; t < keysProfileImage.size(); t++) {
                                    stringToCSVCIVIC.append(objectProfile_image.get(keysProfileImage.get(t))).append(";");
                                }
                            } else {
                                stringToCSVCIVIC.append(objectOrganization.get(keysOrganization.get(d))).append(";");
                            }
                        }
                    } else if (keysUser.get(g).equals("avatars")) {
                        JsonObject objectAvatars = objectUser.get(keysUser.get(g)).getAsJsonObject();
                        List<String> keysAvatars = new ArrayList<>(objectAvatars.keySet());
                        for (int u = 0; u < keysAvatars.size(); u++) {
                            stringToCSVCIVIC.append(objectAvatars.get(keysAvatars.get(u))).append(";");

                        }
                    } else {
                        stringToCSVCIVIC.append(objectUser.get(keysUser.get(g))).append(";");
                    }
                }
            }
        }
    }
}
