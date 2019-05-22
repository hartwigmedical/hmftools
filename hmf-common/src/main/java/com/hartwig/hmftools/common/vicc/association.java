package com.hartwig.hmftools.common.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class association {

    public static StringBuilder readObjectAssociation(@NotNull JsonObject object) {
        //association object
        StringBuilder stringToCSVAssociation = new StringBuilder();

        StringBuilder stringKingdom = new StringBuilder();
        StringBuilder stringDirectParent = new StringBuilder();
        StringBuilder stringClass = new StringBuilder();
        StringBuilder stringSubClass = new StringBuilder();
        StringBuilder stringSuperClass = new StringBuilder();

        StringBuilder stringSource = new StringBuilder();
        StringBuilder stringTerm = new StringBuilder();
        StringBuilder stringDescription = new StringBuilder();
        StringBuilder stringId = new StringBuilder();
        StringBuilder stringUsanStem = new StringBuilder();
        StringBuilder stringApprovedCountries = new StringBuilder();
        StringBuilder stringToxicity = new StringBuilder();

        List<String> keysOfAssocationObject = new ArrayList<>(object.getAsJsonObject("association").keySet());
        // LOGGER.info(keysOfAssocationObject);
        for (int i = 0; i < object.getAsJsonObject("association").keySet().size(); i++) {

            if (keysOfAssocationObject.get(i).equals("drug_labels")) {
                if (i == 0) {
                    JsonElement objectAssociation = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                    stringToCSVAssociation.append(objectAssociation).append(";");
                }
            } else if (keysOfAssocationObject.get(i).equals("description")) {
                if (i == 1 || i == 0) {
                    JsonElement objectAssociation = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                    stringToCSVAssociation.append(objectAssociation).append(";");
                }
            } else if (keysOfAssocationObject.get(i).equals("publication_url")) {
                if (i == 2 || i == 1) {
                    JsonElement objectPublicationUrl = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                    stringToCSVAssociation.append(objectPublicationUrl).append(";");
                } else {
                    stringToCSVAssociation.append(";");
                }
            } else if (keysOfAssocationObject.get(i).equals("source_link")) {
                if (i == 3) {
                    JsonElement objectPublicationUrl = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                    stringToCSVAssociation.append(objectPublicationUrl).append(";");
                } else {
                    stringToCSVAssociation.append(";");
                }
            } else if (keysOfAssocationObject.get(i).equals("variant_name")) {
                if (i == 4 || i == 3 || i == 2) {
                    JsonElement objectPublicationUrl = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                    stringToCSVAssociation.append(objectPublicationUrl).append(";");
                } else {
                    stringToCSVAssociation.append(";");
                }
            } else if (keysOfAssocationObject.get(i).equals("evidence")) {
                if (i == 5 || i == 2 || i == 3 || i == 4 || i == 1) {
                    if (i == 2) {
                        stringToCSVAssociation.append(";;;");
                    }
                    JsonArray arrayEvidence = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i)).getAsJsonArray();
                    for (int x = 0; x < arrayEvidence.size(); x++) {
                        JsonObject objectEvidence = (JsonObject) arrayEvidence.get(x);
                        List<String> keysEvidence = new ArrayList<>(objectEvidence.keySet());
                        for (int d = 0; d < keysEvidence.size(); d++) {
                            if (keysEvidence.get(d).equals("info")) {
                                if (objectEvidence.get(keysEvidence.get(d)).isJsonNull()) {
                                    stringToCSVAssociation.append(";;");
                                } else {
                                    JsonObject objectInfo = objectEvidence.get(keysEvidence.get(d)).getAsJsonObject();

                                    List<String> keysInfo = new ArrayList<>(objectInfo.keySet());
                                    for (int z = 0; z < keysInfo.size(); z++) {
                                        stringToCSVAssociation.append(objectInfo.get(keysInfo.get(z))).append(";");
                                    }
                                }
                            } else if (keysEvidence.get(d).equals("evidenceType")) {
                                JsonObject objectEvidenceType = objectEvidence.get(keysEvidence.get(d)).getAsJsonObject();
                                List<String> keysEvidenceType = new ArrayList<>(objectEvidenceType.keySet());
                                for (int z = 0; z < keysEvidenceType.size(); z++) {
                                    stringToCSVAssociation.append(objectEvidenceType.get(keysEvidenceType.get(z))).append(";");
                                    if (i == 2 || i == 4) {
                                        stringToCSVAssociation.append(";");
                                    }
                                }
                            } else if (keysEvidence.get(d).equals("description")) {
                                stringToCSVAssociation.append(objectEvidence.get(keysEvidence.get(d))).append(";");
                            }
                        }
                    }
                }
            } else if (keysOfAssocationObject.get(i).equals("environmentalContexts")) {
                if (i == 6 || i == 3 || i == 4 || i == 2) {
                    JsonArray arrayEvidence = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i)).getAsJsonArray();
                    for (int u = 0; u < arrayEvidence.size(); u++) {
                        JsonObject objectEnvironmentalContexts = (JsonObject) arrayEvidence.get(u);
                        List<String> keysEnvironmentalContexts = new ArrayList<>(objectEnvironmentalContexts.keySet());
                        for (int g = 0; g < keysEnvironmentalContexts.size(); g++) {
                            if (keysEnvironmentalContexts.get(g).equals("taxonomy")) {
                                JsonObject objectTaxonomy = objectEnvironmentalContexts.getAsJsonObject(keysEnvironmentalContexts.get(g));
                                List<String> keysTaxonomy = new ArrayList<>(objectTaxonomy.keySet());
                                for (int p = 0; p < keysTaxonomy.size(); p++) {
                                    if (keysTaxonomy.get(p).equals("kingdom")) {
                                        JsonElement kingdom = objectTaxonomy.get(keysTaxonomy.get(p));
                                        stringKingdom.append(kingdom).append(",");
                                    } else if (keysTaxonomy.get(p).equals("direct-parent")) {
                                        JsonElement directParent = objectTaxonomy.get(keysTaxonomy.get(p));
                                        stringDirectParent.append(directParent).append(",");
                                    } else if (keysTaxonomy.get(p).equals("class")) {
                                        JsonElement classTaxanomy = objectTaxonomy.get(keysTaxonomy.get(p));
                                        stringClass.append(classTaxanomy).append(",");
                                    } else if (keysTaxonomy.get(p).equals("subclass")) {
                                        JsonElement subClass = objectTaxonomy.get(keysTaxonomy.get(p));
                                        stringSubClass.append(subClass).append(",");
                                    } else if (keysTaxonomy.get(p).equals("superclass")) {
                                        JsonElement superClass = objectTaxonomy.get(keysTaxonomy.get(p));
                                        stringSuperClass.append(superClass).append(",");
                                    }
                                }
                            } else {
                                if (keysEnvironmentalContexts.get(g).equals("source")) {
                                    JsonElement source = objectEnvironmentalContexts.get(keysEnvironmentalContexts.get(g));
                                    stringSource.append(source).append(",");
                                } else if (keysEnvironmentalContexts.get(g).equals("term")) {
                                    JsonElement term = objectEnvironmentalContexts.get(keysEnvironmentalContexts.get(g));
                                    stringTerm.append(term).append(",");
                                } else if (keysEnvironmentalContexts.get(g).equals("description")) {
                                    JsonElement description = objectEnvironmentalContexts.get(keysEnvironmentalContexts.get(g));
                                    stringDescription.append(description).append(",");
                                } else if (keysEnvironmentalContexts.get(g).equals("id")) {
                                    JsonElement id = objectEnvironmentalContexts.get(keysEnvironmentalContexts.get(g));
                                    stringId.append(id).append(",");
                                } else if (keysEnvironmentalContexts.get(g).equals("usan_stem")) {
                                    JsonElement usanStem = objectEnvironmentalContexts.get(keysEnvironmentalContexts.get(g));
                                    stringUsanStem.append(usanStem).append(",");
                                } else if (keysEnvironmentalContexts.get(g).equals("approved_countries")) {
                                    JsonElement approvedCountries = objectEnvironmentalContexts.get(keysEnvironmentalContexts.get(g));
                                    stringApprovedCountries.append(approvedCountries);
                                } else if (keysEnvironmentalContexts.get(g).equals("toxicity")) {
                                    JsonElement toxity = objectEnvironmentalContexts.get(keysEnvironmentalContexts.get(g));
                                    stringToxicity.append(toxity);
                                }
                            }
                        }
                    }
                }
                stringToCSVAssociation.append(stringTerm)
                        .append(";")
                        .append(stringDescription)
                        .append(";")
                        .append(stringKingdom)
                        .append(";")
                        .append(stringDirectParent)
                        .append(";")
                        .append(stringClass)
                        .append(";")
                        .append(stringSubClass)
                        .append(";")
                        .append(stringSuperClass)
                        .append(";")
                        .append(stringSource)
                        .append(";")
                        .append(stringUsanStem)
                        .append(";")
                        .append(stringToxicity)
                        .append(";")
                        .append(stringApprovedCountries)
                        .append(";")
                        .append(stringId)
                        .append(";");
            } else if (keysOfAssocationObject.get(i).equals("evidence_label")) {
                if (i == 7 || i == 4 || i == 5 || i == 3) {
                    JsonElement objectAssociation = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                    stringToCSVAssociation.append(objectAssociation).append(";");
                }
            } else if (keysOfAssocationObject.get(i).equals("phenotype")) {
                if (i == 8 || i == 5 || i == 6 || i == 4) {
                    JsonObject objectPhenotype = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i)).getAsJsonObject();
                    List<String> keysPhenotype = new ArrayList<>(objectPhenotype.keySet());
                    for (int y = 0; y < keysPhenotype.size(); y++) {
                        if (keysPhenotype.get(y).equals("type")) {
                            JsonObject objectInfo = objectPhenotype.get(keysPhenotype.get(y)).getAsJsonObject();
                            List<String> keysInfo = new ArrayList<>(objectInfo.keySet());
                            for (int z = 0; z < keysInfo.size(); z++) {
                                stringToCSVAssociation.append(objectInfo.get(keysInfo.get(z))).append(";");
                            }
                        } else {
                            stringToCSVAssociation.append(objectPhenotype.get(keysPhenotype.get(y))).append(";");
                            if (i == 5 || i == 7 && keysPhenotype.get(y).equals("family")) {
                                stringToCSVAssociation.append(";");
                            }
                        }
                    }
                }
            } else if (keysOfAssocationObject.get(i).equals("evidence_level")) {
                if (i == 9 || i == 6 || i == 7) {
                    JsonElement objectAssociation = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                    stringToCSVAssociation.append(objectAssociation).append(";");
                }
            } else if (keysOfAssocationObject.get(i).equals("response_type")) {
                if (i == 10 || i == 8 || i == 7) {
                    JsonElement objectAssociation = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                    stringToCSVAssociation.append(objectAssociation).append(";");
                }
            } else if (keysOfAssocationObject.get(i).equals("oncogenic")) {
                if (i == 5 || i == 7) {
                    JsonElement objectAssociation = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                    stringToCSVAssociation.append(objectAssociation).append(";");
                }
            }
        }
        return stringToCSVAssociation;
    }
}
