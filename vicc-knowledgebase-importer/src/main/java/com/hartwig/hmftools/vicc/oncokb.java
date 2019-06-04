package com.hartwig.hmftools.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class oncokb {

    public static StringBuilder readObjectOncokb(@NotNull JsonObject object) {
        //ONCOKB object
        StringBuilder stringToCSVOncoKb = new StringBuilder();
        StringBuilder info = new StringBuilder();
        StringBuilder extra = new StringBuilder();
        StringBuilder biologicalInfo = new StringBuilder();
        StringBuilder clinicallInfo = new StringBuilder();
        StringBuilder stringLink = new StringBuilder();
        StringBuilder stringText = new StringBuilder();

        if (object.getAsJsonObject("oncokb") != null) {
            List<String> keysOfoncoKb = new ArrayList<>(object.getAsJsonObject("oncokb").keySet());
            for (int x = 0; x < keysOfoncoKb.size(); x++) {
                JsonObject pmkbObject = object.getAsJsonObject("oncokb").get(keysOfoncoKb.get(x)).getAsJsonObject();
                List<String> keysOfBiological = new ArrayList<>(pmkbObject.keySet());
                for (int y = 0; y < keysOfBiological.size(); y++) {
                    if (keysOfoncoKb.get(x).equals("biological")) {
                        if (keysOfBiological.get(y).equals("mutationEffectPmids")) {
                            info.append(pmkbObject.get(keysOfBiological.get(y))).append(";");
                        } else if (keysOfBiological.get(y).equals("Isoform")) {
                            info.append(pmkbObject.get(keysOfBiological.get(y))).append(";");
                            info.append(";;;");
                        } else if (keysOfBiological.get(y).equals("mutationEffectAbstracts")) {
                            biologicalInfo.append(pmkbObject.get(keysOfBiological.get(y))).append(";");
                            biologicalInfo.append(";;;;;;;;");
                        } else if (keysOfBiological.get(y).equals("Entrez Gene ID") || keysOfBiological.get(y).equals("oncogenic")
                                || keysOfBiological.get(y).equals("mutationEffect") || keysOfBiological.get(y).equals("RefSeq")
                                || keysOfBiological.get(y).equals("gene")) {
                            biologicalInfo.append(pmkbObject.get(keysOfBiological.get(y))).append(";");
                        }

                    } else if (keysOfoncoKb.get(x).equals("clinical")) {
                        if (keysOfBiological.get(y).equals("RefSeq")) {
                            info.append(";;");
                            info.append(pmkbObject.get(keysOfBiological.get(y))).append(";");
                        } else if (keysOfBiological.get(y).equals("level") || keysOfBiological.get(y).equals("Isoform")) {
                            info.append(pmkbObject.get(keysOfBiological.get(y))).append(";");
                        } else if (keysOfBiological.get(y).equals("drugAbstracts")) {
                            JsonArray drugs = pmkbObject.get(keysOfBiological.get(y)).getAsJsonArray();
                            for (int v = 0; v < drugs.size(); v++) {
                                JsonObject objectDrugs = (JsonObject) drugs.get(v);
                                List<String> keysDrugs = new ArrayList<>(objectDrugs.keySet());
                                for (int u = 0; u < keysDrugs.size(); u++) {
                                    if (keysDrugs.get(u).equals("text")) {
                                        JsonElement text = objectDrugs.get(keysDrugs.get(u));
                                        stringText.append(text).append(",");
                                    } else if (keysDrugs.get(u).equals("link")) {
                                        JsonElement link = objectDrugs.get(keysDrugs.get(u));
                                        stringLink.append(link).append(",");
                                    }
                                }
                            }
                            clinicallInfo.append(stringText).append(";").append(stringLink).append(";");
                        } else if (keysOfBiological.get(y).equals("Entrez Gene ID")) {
                            clinicallInfo.append(";;;;;;");
                            clinicallInfo.append(pmkbObject.get(keysOfBiological.get(y))).append(";");
                        } else if (keysOfBiological.get(y).equals("drugPmids") || keysOfBiological.get(y).equals("cancerType")
                                || keysOfBiological.get(y).equals("drug") || keysOfBiological.get(y).equals("gene") || keysOfBiological.get(
                                y).equals("level_label")) {
                            clinicallInfo.append(pmkbObject.get(keysOfBiological.get(y))).append(";");
                        }
                    }
                    if (keysOfBiological.get(y).equals("variant")) {
                        JsonObject oncokbVariant = pmkbObject.get(keysOfBiological.get(y)).getAsJsonObject();
                        List<String> keysOfVariant = new ArrayList<>(oncokbVariant.keySet());
                        for (int z = 0; z < keysOfVariant.size(); z++) {
                            if (keysOfVariant.get(z).equals("consequence")) {
                                JsonObject oncokbConsequence = oncokbVariant.get(keysOfVariant.get(z)).getAsJsonObject();
                                List<String> keysOfConsequence = new ArrayList<>(oncokbConsequence.keySet());
                                for (int v = 0; v < keysOfConsequence.size(); v++) {
                                    extra.append(oncokbConsequence.get(keysOfConsequence.get(v))).append(";");
                                }
                            } else if (keysOfVariant.get(z).equals("gene")) {
                                JsonObject oncokbGene = oncokbVariant.get(keysOfVariant.get(z)).getAsJsonObject();
                                List<String> keysOfGene = new ArrayList<>(oncokbGene.keySet());
                                for (int i = 0; i < keysOfGene.size(); i++) {
                                    extra.append(oncokbGene.get(keysOfGene.get(i))).append(";");
                                }
                            } else {
                                extra.append(oncokbVariant.get(keysOfVariant.get(z))).append(";");
                            }
                        }
                    } else if (!keysOfBiological.get(y).equals("mutationEffectPmids") && !keysOfBiological.get(y).equals("Isoform")
                            && !keysOfBiological.get(y).equals("RefSeq") && !keysOfBiological.get(y).equals("level")
                            && !keysOfBiological.get(y).equals("Entrez Gene ID") && !keysOfBiological.get(y).equals("drugPmids")
                            && !keysOfBiological.get(y).equals("cancerType") && !keysOfBiological.get(y).equals("drug") && !keysOfBiological
                            .get(y)
                            .equals("gene") && !keysOfBiological.get(y).equals("level_label") && !keysOfBiological.get(y)
                            .equals("oncogenic") && !keysOfBiological.get(y).equals("mutationEffect") && !keysOfBiological.get(y)
                            .equals("mutationEffectAbstracts") && !keysOfBiological.get(y).equals("drugAbstracts")) {
                        extra.append(pmkbObject.get(keysOfBiological.get(y))).append(";");
                    }
                }
            }
            stringToCSVOncoKb.append(info).append(extra).append(biologicalInfo).append(clinicallInfo);
        } else {
            stringToCSVOncoKb.append(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;");
        }
        return stringToCSVOncoKb;
    }
}
