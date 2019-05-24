package com.hartwig.hmftools.common.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class molecularMatch {

    public static StringBuilder readObjectMolecularMatch(@NotNull JsonObject object) {
        //MolecularMatch object
        StringBuilder stringToCSVMolecularMatch = new StringBuilder();

        StringBuilder stringPriority = new StringBuilder();
        StringBuilder stringCompositeKey = new StringBuilder();
        StringBuilder stringSuppress = new StringBuilder();
        StringBuilder stringFilterType = new StringBuilder();
        StringBuilder stringTerm = new StringBuilder();
        StringBuilder stringPrimary = new StringBuilder();
        StringBuilder stringFacet = new StringBuilder();
        StringBuilder stringValid = new StringBuilder();
        StringBuilder stringCustom = new StringBuilder();

        StringBuilder stringCount = new StringBuilder();
        StringBuilder stringPercent = new StringBuilder();
        StringBuilder stringStudyId = new StringBuilder();
        StringBuilder stringSamples = new StringBuilder();
        StringBuilder stringMolecular = new StringBuilder();
        StringBuilder stringCondition = new StringBuilder();

        StringBuilder stringAminoAcidChange = new StringBuilder();
        StringBuilder stringCompositeKey1 = new StringBuilder();
        StringBuilder stringIntronNumber = new StringBuilder();
        StringBuilder stringRef = new StringBuilder();
        StringBuilder stringExonNumber = new StringBuilder();
        StringBuilder stringSupress = new StringBuilder();
        StringBuilder stringStop = new StringBuilder();
        StringBuilder stringCustom1 = new StringBuilder();
        StringBuilder stringStart = new StringBuilder();
        StringBuilder stringChr = new StringBuilder();
        StringBuilder stringStrand = new StringBuilder();
        StringBuilder stringAlt = new StringBuilder();
        StringBuilder stringValidated = new StringBuilder();
        StringBuilder stringTranscript = new StringBuilder();
        StringBuilder stringCdna = new StringBuilder();
        StringBuilder stringReferenceGenome = new StringBuilder();

        StringBuilder stringAminoAcidChangeRefLocation = new StringBuilder();
        StringBuilder stringTxSitesRefLocation = new StringBuilder();
        StringBuilder stringExonNumberRefLocation = new StringBuilder();
        StringBuilder stringIntronNumberRefLocation = new StringBuilder();
        StringBuilder stringTranscriptRefLocation = new StringBuilder();
        StringBuilder stringCdnaRefLocation = new StringBuilder();

        StringBuilder stringNameSources = new StringBuilder();
        StringBuilder stringSuppressSources = new StringBuilder();
        StringBuilder stringpubIdSources = new StringBuilder();
        StringBuilder stringSubTypeSources = new StringBuilder();
        StringBuilder stringValidSources = new StringBuilder();
        StringBuilder stringLinkSources = new StringBuilder();
        StringBuilder stringYearSources = new StringBuilder();
        StringBuilder stringTypeSources = new StringBuilder();
        StringBuilder stringIdSources = new StringBuilder();

        StringBuilder stringTier = new StringBuilder();
        StringBuilder stringStep = new StringBuilder();
        StringBuilder stringMessage = new StringBuilder();
        StringBuilder stringSuccess = new StringBuilder();

        StringBuilder stringPriorityTags = new StringBuilder();
        StringBuilder stringCompositeKeyTags = new StringBuilder();
        StringBuilder stringSuppressTags = new StringBuilder();
        StringBuilder stringFilterTypeTags = new StringBuilder();
        StringBuilder stringTermTags = new StringBuilder();
        StringBuilder stringPrimaryTags = new StringBuilder();
        StringBuilder stringFacetTags = new StringBuilder();
        StringBuilder stringValidTags = new StringBuilder();
        StringBuilder stringCustomTags = new StringBuilder();

        if (object.getAsJsonObject("molecularmatch") != null) {
            List<String> keysOfMolecularMatch = new ArrayList<>(object.getAsJsonObject("molecularmatch").keySet());
            for (int x = 0; x < keysOfMolecularMatch.size(); x++) {
                if (keysOfMolecularMatch.get(x).equals("criteriaUnmet")) {
                    JsonArray molecluarMatchArray =
                            object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchArray.size(); i++) {
                        JsonObject objectVriteriaUnmet = (JsonObject) molecluarMatchArray.get(i);
                        List<String> keysCriteriaUnmet = new ArrayList<>(objectVriteriaUnmet.keySet());
                        for (int y = 0; y < keysCriteriaUnmet.size(); y++) {
                            if (keysCriteriaUnmet.get(y).equals("priority")) {
                                JsonElement priority = objectVriteriaUnmet.get(keysCriteriaUnmet.get(y));
                                stringPriority.append(priority).append(",");
                            } else if (keysCriteriaUnmet.get(y).equals("compositeKey")) {
                                JsonElement compositeKey = objectVriteriaUnmet.get(keysCriteriaUnmet.get(y));
                                stringCompositeKey.append(compositeKey).append(",");
                            } else if (keysCriteriaUnmet.get(y).equals("suppress")) {
                                JsonElement suppress = objectVriteriaUnmet.get(keysCriteriaUnmet.get(y));
                                stringSuppress.append(suppress).append(",");
                            } else if (keysCriteriaUnmet.get(y).equals("filterType")) {
                                JsonElement filterType = objectVriteriaUnmet.get(keysCriteriaUnmet.get(y));
                                stringFilterType.append(filterType).append(",");
                            } else if (keysCriteriaUnmet.get(y).equals("term")) {
                                JsonElement term = objectVriteriaUnmet.get(keysCriteriaUnmet.get(y));
                                stringTerm.append(term).append(",");
                            } else if (keysCriteriaUnmet.get(y).equals("primary")) {
                                JsonElement primary = objectVriteriaUnmet.get(keysCriteriaUnmet.get(y));
                                stringPrimary.append(primary).append(",");
                            } else if (keysCriteriaUnmet.get(y).equals("facet")) {
                                JsonElement facet = objectVriteriaUnmet.get(keysCriteriaUnmet.get(y));
                                stringFacet.append(facet).append(",");
                            } else if (keysCriteriaUnmet.get(y).equals("valid")) {
                                JsonElement valid = objectVriteriaUnmet.get(keysCriteriaUnmet.get(y));
                                stringValid.append(valid).append(",");
                            } else if (keysCriteriaUnmet.get(y).equals("custom")) {
                                JsonElement custom = objectVriteriaUnmet.get(keysCriteriaUnmet.get(y));
                                stringCustom.append(custom).append(",");
                            }
                        }
                    }
                    stringToCSVMolecularMatch.append(stringPriority)
                            .append(";")
                            .append(stringCompositeKey)
                            .append(";")
                            .append(stringSuppress)
                            .append(";")
                            .append(stringFilterType)
                            .append(";")
                            .append(stringTerm)
                            .append(";")
                            .append(stringPrimary)
                            .append(";")
                            .append(stringFacet)
                            .append(";")
                            .append(stringValid)
                            .append(";")
                            .append(stringCustom)
                            .append(";");
                } else if (keysOfMolecularMatch.get(x).equals("prevalence")) {
                    JsonArray molecluarMatchArray =
                            object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchArray.size(); i++) {
                        JsonObject objectPrevalence = (JsonObject) molecluarMatchArray.get(i);
                        List<String> keysPrevalence = new ArrayList<>(objectPrevalence.keySet());
                        for (int y = 0; y < keysPrevalence.size(); y++) {
                            if (keysPrevalence.get(y).equals("count")) {
                                JsonElement count = objectPrevalence.get(keysPrevalence.get(y));
                                stringCount.append(count).append(",");
                            } else if (keysPrevalence.get(y).equals("percent")) {
                                JsonElement percent = objectPrevalence.get(keysPrevalence.get(y));
                                stringPercent.append(percent).append(",");
                            } else if (keysPrevalence.get(y).equals("studyId")) {
                                JsonElement studyId = objectPrevalence.get(keysPrevalence.get(y));
                                stringStudyId.append(studyId).append(",");
                            } else if (keysPrevalence.get(y).equals("samples")) {
                                JsonElement samples = objectPrevalence.get(keysPrevalence.get(y));
                                stringSamples.append(samples).append(",");
                            } else if (keysPrevalence.get(y).equals("molecular")) {
                                JsonElement molecular = objectPrevalence.get(keysPrevalence.get(y));
                                stringMolecular.append(molecular).append(",");
                            } else if (keysPrevalence.get(y).equals("condition")) {
                                JsonElement condition = objectPrevalence.get(keysPrevalence.get(y));
                                stringCondition.append(condition).append(",");
                            }
                        }
                    }
                    stringToCSVMolecularMatch.append(stringCount)
                            .append(";")
                            .append(stringPercent)
                            .append(";")
                            .append(stringStudyId)
                            .append(";")
                            .append(stringSamples)
                            .append(";")
                            .append(stringMolecular)
                            .append(";")
                            .append(stringCondition)
                            .append(";");
                } else if (keysOfMolecularMatch.get(x).equals("mutations")) {
                    JsonArray molecluarMatchArray =
                            object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchArray.size(); i++) {
                        JsonObject objectMutations = (JsonObject) molecluarMatchArray.get(i);
                        List<String> keysMutations = new ArrayList<>(objectMutations.keySet());
                        for (int u = 0; u < keysMutations.size(); u++) {
                            if (keysMutations.get(u).equals("transcriptConsequence")) {
                                JsonArray arrayTranscriptConsequence = objectMutations.get("transcriptConsequence").getAsJsonArray();
                                for (int g = 0; g < arrayTranscriptConsequence.size(); g++) {
                                    JsonObject objectTranscriptConsequence = (JsonObject) arrayTranscriptConsequence.get(g);
                                    List<String> keysTranscriptConsequence = new ArrayList<>(objectTranscriptConsequence.keySet());
                                    for (int p = 0; p < keysTranscriptConsequence.size(); p++) {
                                        if (keysTranscriptConsequence.get(p).equals("amino_acid_change")) {
                                            JsonElement aminoAcidChange = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringAminoAcidChange.append(aminoAcidChange).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("compositeKey")) {
                                            JsonElement compositeKey = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringCompositeKey1.append(compositeKey).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("intronNumber")) {
                                            JsonElement intronNumber = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringIntronNumber.append(intronNumber).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("ref")) {
                                            JsonElement ref = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringRef.append(ref).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("exonNumber")) {
                                            JsonElement exonNumber = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringExonNumber.append(exonNumber).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("suppress")) {
                                            JsonElement suppress = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringSupress.append(suppress).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("stop")) {
                                            JsonElement stop = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringStop.append(stop).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("custom")) {
                                            JsonElement custom = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringCustom1.append(custom).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("start")) {
                                            JsonElement start = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringStart.append(start).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("chr")) {
                                            JsonElement chr = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringChr.append(chr).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("strand")) {
                                            JsonElement strand = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringStrand.append(strand).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("alt")) {
                                            JsonElement alt = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringAlt.append(alt).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("validated")) {
                                            JsonElement validated = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringValidated.append(validated).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("transcript")) {
                                            JsonElement transcript = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringTranscript.append(transcript).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("cdna")) {
                                            JsonElement cdna = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringCdna.append(cdna).append(",");
                                        } else if (keysTranscriptConsequence.get(p).equals("referenceGenome")) {
                                            JsonElement referenceGenome = objectTranscriptConsequence.get(keysTranscriptConsequence.get(p));
                                            stringReferenceGenome.append(referenceGenome).append(",");
                                        }
                                    }
                                }
                                stringToCSVMolecularMatch.append(stringAminoAcidChange)
                                        .append(";")
                                        .append(stringCompositeKey1)
                                        .append(";")
                                        .append(stringIntronNumber)
                                        .append(";")
                                        .append(stringRef)
                                        .append(";")
                                        .append(stringExonNumber)
                                        .append(";")
                                        .append(stringSupress)
                                        .append(";")
                                        .append(stringStop)
                                        .append(";")
                                        .append(stringCustom1)
                                        .append(";")
                                        .append(stringStart)
                                        .append(";")
                                        .append(stringChr)
                                        .append(";")
                                        .append(stringStrand)
                                        .append(";")
                                        .append(stringAlt)
                                        .append(";")
                                        .append(stringValidated)
                                        .append(";")
                                        .append(stringTranscript)
                                        .append(";")
                                        .append(stringCdna)
                                        .append(";")
                                        .append(stringReferenceGenome)
                                        .append(";");
                            } else if (keysMutations.get(u).equals("GRCh37_location")) {
                                JsonArray arrayLocationRefGenome = objectMutations.get("GRCh37_location").getAsJsonArray();
                                for (int r = 0; r < arrayLocationRefGenome.size(); r++) {
                                    JsonObject objectLocationRefGenome = (JsonObject) arrayLocationRefGenome.get(r);
                                    List<String> keysLocationRefGenome = new ArrayList<>(objectLocationRefGenome.keySet());
                                    for (int g = 0; g < keysLocationRefGenome.size(); g++) {
                                        if (keysLocationRefGenome.get(g).equals("transcript_consequences")) {
                                            JsonArray arrayTranscriptConsequence =
                                                    objectLocationRefGenome.get("transcript_consequences").getAsJsonArray();
                                            for (int e = 0; e < arrayTranscriptConsequence.size(); e++) {
                                                JsonObject objectTranscriptConsequence = (JsonObject) arrayTranscriptConsequence.get(e);
                                                List<String> keysTranscriptConsequence =
                                                        new ArrayList<>(objectTranscriptConsequence.keySet());
                                                for (int o = 0; o < keysTranscriptConsequence.size(); o++) {
                                                    if (keysTranscriptConsequence.get(o).equals("amino_acid_change")) {
                                                        JsonElement aminoAcidChange =
                                                                objectTranscriptConsequence.get(keysTranscriptConsequence.get(o));
                                                        stringAminoAcidChangeRefLocation.append(aminoAcidChange).append(",");
                                                    } else if (keysTranscriptConsequence.get(o).equals("txSites")) {
                                                        JsonElement txSites =
                                                                objectTranscriptConsequence.get(keysTranscriptConsequence.get(o));
                                                        stringTxSitesRefLocation.append(txSites).append(",");
                                                    } else if (keysTranscriptConsequence.get(o).equals("exonNumber")) {
                                                        JsonElement exonNumber =
                                                                objectTranscriptConsequence.get(keysTranscriptConsequence.get(o));
                                                        stringExonNumberRefLocation.append(exonNumber).append(",");
                                                    } else if (keysTranscriptConsequence.get(o).equals("intronNumber")) {
                                                        JsonElement intronNumber =
                                                                objectTranscriptConsequence.get(keysTranscriptConsequence.get(o));
                                                        stringIntronNumberRefLocation.append(intronNumber).append(",");
                                                    } else if (keysTranscriptConsequence.get(o).equals("transcript")) {
                                                        JsonElement transcript =
                                                                objectTranscriptConsequence.get(keysTranscriptConsequence.get(o));
                                                        stringTranscriptRefLocation.append(transcript).append(",");
                                                    } else if (keysTranscriptConsequence.get(o).equals("cdna")) {
                                                        JsonElement cdna =
                                                                objectTranscriptConsequence.get(keysTranscriptConsequence.get(o));
                                                        stringCdnaRefLocation.append(cdna).append(",");
                                                    }
                                                }
                                            }
                                            stringToCSVMolecularMatch.append(stringAminoAcidChangeRefLocation)
                                                    .append(";")
                                                    .append(stringTxSitesRefLocation)
                                                    .append(";")
                                                    .append(stringExonNumberRefLocation)
                                                    .append(";")
                                                    .append(stringIntronNumberRefLocation)
                                                    .append(";")
                                                    .append(stringTranscriptRefLocation)
                                                    .append(";")
                                                    .append(stringCdnaRefLocation)
                                                    .append(";");
                                        } else {
                                            JsonElement elementChRCHLocation = objectLocationRefGenome.get(keysLocationRefGenome.get(g));
                                            stringToCSVMolecularMatch.append(elementChRCHLocation).append(";");
                                        }
                                    }
                                }
                            } else {
                                JsonArray arrays = objectMutations.get("sources").getAsJsonArray();
                                stringToCSVMolecularMatch.append(arrays).append(";");
                            }
                        }
                    }
                } else if (keysOfMolecularMatch.get(x).equals("sources")) {
                    JsonArray molecluarMatchArray =
                            object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchArray.size(); i++) {
                        JsonObject objectSources = (JsonObject) molecluarMatchArray.get(i);
                        List<String> keysSources = new ArrayList<>(objectSources.keySet());
                        for (int u = 0; u < keysSources.size(); u++) {
                            if (keysSources.get(u).equals("name")) {
                                JsonElement name = objectSources.get(keysSources.get(u));
                                stringNameSources.append(name).append(",");
                            } else if (keysSources.get(u).equals("suppress")) {
                                JsonElement suppress = objectSources.get(keysSources.get(u));
                                stringSuppressSources.append(suppress).append(",");
                            } else if (keysSources.get(u).equals("pubId")) {
                                JsonElement pubId = objectSources.get(keysSources.get(u));
                                stringpubIdSources.append(pubId).append(",");
                            } else if (keysSources.get(u).equals("subType")) {
                                JsonElement subType = objectSources.get(keysSources.get(u));
                                stringSubTypeSources.append(subType).append(",");
                            } else if (keysSources.get(u).equals("valid")) {
                                JsonElement valid = objectSources.get(keysSources.get(u));
                                stringValidSources.append(valid).append(",");
                            } else if (keysSources.get(u).equals("link")) {
                                JsonElement link = objectSources.get(keysSources.get(u));
                                stringLinkSources.append(link).append(",");
                            } else if (keysSources.get(u).equals("year")) {
                                JsonElement year = objectSources.get(keysSources.get(u));
                                stringYearSources.append(year).append(",");
                            } else if (keysSources.get(u).equals("type")) {
                                JsonElement type = objectSources.get(keysSources.get(u));
                                stringTypeSources.append(type).append(",");
                            } else if (keysSources.get(u).equals("id")) {
                                JsonElement id = objectSources.get(keysSources.get(u));
                                stringIdSources.append(id).append(",");
                            }

                        }
                    }
                    stringToCSVMolecularMatch.append(stringNameSources)
                            .append(";")
                            .append(stringSuppressSources)
                            .append(";")
                            .append(stringpubIdSources)
                            .append(";")
                            .append(stringSubTypeSources)
                            .append(";")
                            .append(stringValidSources)
                            .append(";")
                            .append(stringLinkSources)
                            .append(";")
                            .append(stringYearSources)
                            .append(";")
                            .append(stringTypeSources)
                            .append(";")
                            .append(stringIdSources)
                            .append(";");
                } else if (keysOfMolecularMatch.get(x).equals("ast")) {
                    JsonObject molecularMatchObject =
                            object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)).getAsJsonObject();
                    List<String> keysast = new ArrayList<>(molecularMatchObject.keySet());
                    for (int q = 0; q < keysast.size(); q++) {
                        if (keysast.get(q).equals("right")) {
                            List<String> keysright = new ArrayList<>(molecularMatchObject.get(keysast.get(q)).getAsJsonObject().keySet());
                            for (int f = 0; f < keysright.size(); f++) {
                                stringToCSVMolecularMatch.append(keysright).append(";");
                            }
                        } else if (keysast.get(q).equals("left")) {
                            List<String> keysleft = new ArrayList<>(molecularMatchObject.get(keysast.get(q)).getAsJsonObject().keySet());
                            for (int f = 0; f < keysleft.size(); f++) {
                                stringToCSVMolecularMatch.append(keysleft).append(";");
                            }
                        } else {
                            stringToCSVMolecularMatch.append(keysast).append(";");
                        }
                    }
                } else if (keysOfMolecularMatch.get(x).equals("variantInfo")) {
                    JsonArray molecluarMatchArray =
                            object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchArray.size(); i++) {
                        JsonObject objectVariantInfo = (JsonObject) molecluarMatchArray.get(i);
                        List<String> keysVariantInfo = new ArrayList<>(objectVariantInfo.keySet());
                        for (int u = 0; u < keysVariantInfo.size(); u++) {
                            if (keysVariantInfo.get(u).equals("consequences")) {
                                stringToCSVMolecularMatch.append(objectVariantInfo.get(keysVariantInfo.get(u))).append(";");
                            } else if (keysVariantInfo.get(u).equals("locations")) {
                                JsonArray arrayLocations = objectVariantInfo.get(keysVariantInfo.get(u)).getAsJsonArray();
                                for (int h = 0; h < arrayLocations.size(); h++) {
                                    JsonObject objectLocations = (JsonObject) arrayLocations.get(h);
                                    List<String> keysLocations = new ArrayList<>(objectLocations.keySet());
                                    for (int d = 0; d < keysLocations.size(); d++) {
                                        stringToCSVMolecularMatch.append(objectLocations.get(keysLocations.get(d))).append(";");
                                    }
                                }
                            }
                        }
                    }
                } else if (keysOfMolecularMatch.get(x).equals("tierExplanation")) {
                    JsonArray molecluarMatchArray =
                            object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchArray.size(); i++) {
                        JsonObject objectTierexplanation = (JsonObject) molecluarMatchArray.get(i);
                        List<String> keysTierExplanation = new ArrayList<>(objectTierexplanation.keySet());
                        for (int s = 0; s < keysTierExplanation.size(); s++) {
                            if (keysTierExplanation.get(s).equals("tier")) {
                                JsonElement tier = objectTierexplanation.get(keysTierExplanation.get(s));
                                stringTier.append(tier).append(",");
                            } else if (keysTierExplanation.get(s).equals("step")) {
                                JsonElement step = objectTierexplanation.get(keysTierExplanation.get(s));
                                stringStep.append(step).append(",");
                            } else if (keysTierExplanation.get(s).equals("message")) {
                                JsonElement message = objectTierexplanation.get(keysTierExplanation.get(s));
                                stringMessage.append(message).append(",");
                            } else if (keysTierExplanation.get(s).equals("success")) {
                                JsonElement success = objectTierexplanation.get(keysTierExplanation.get(s));
                                stringSuccess.append(success).append(",");
                            }
                        }
                    }
                    stringToCSVMolecularMatch.append(stringTier)
                            .append(";")
                            .append(stringStep)
                            .append(";")
                            .append(stringMessage)
                            .append(";")
                            .append(stringSuccess)
                            .append(";");
                } else if (keysOfMolecularMatch.get(x).equals("tags")) {
                    JsonArray molecluarMatchArray =
                            object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchArray.size(); i++) {
                        JsonObject objectTags = (JsonObject) molecluarMatchArray.get(i);
                        List<String> keysTags = new ArrayList<>(objectTags.keySet());
                        for (int s = 0; s < keysTags.size(); s++) {
                            if (keysTags.get(s).equals("priority")) {
                                JsonElement priority = objectTags.get(keysTags.get(s));
                                stringPriorityTags.append(priority).append(",");
                            } else if (keysTags.get(s).equals("compositeKey")) {
                                JsonElement compositeKey = objectTags.get(keysTags.get(s));
                                stringCompositeKeyTags.append(compositeKey).append(",");
                            } else if (keysTags.get(s).equals("suppress")) {
                                JsonElement suppress = objectTags.get(keysTags.get(s));
                                stringSuppressTags.append(suppress).append(",");
                            } else if (keysTags.get(s).equals("filterType")) {
                                JsonElement filterType = objectTags.get(keysTags.get(s));
                                stringFilterTypeTags.append(filterType).append(",");
                            } else if (keysTags.get(s).equals("term")) {
                                JsonElement term = objectTags.get(keysTags.get(s));
                                stringTermTags.append(term).append(",");
                            } else if (keysTags.get(s).equals("primary")) {
                                JsonElement primary = objectTags.get(keysTags.get(s));
                                stringPrimaryTags.append(primary).append(",");
                            } else if (keysTags.get(s).equals("facet")) {
                                JsonElement facet = objectTags.get(keysTags.get(s));
                                stringFacetTags.append(facet).append(",");
                            } else if (keysTags.get(s).equals("valid")) {
                                JsonElement valid = objectTags.get(keysTags.get(s));
                                stringValidTags.append(valid).append(",");
                            } else if (keysTags.get(s).equals("custom")) {
                                JsonElement custom = objectTags.get(keysTags.get(s));
                                stringCustomTags.append(custom).append(",");
                            }
                        }
                    }
                    stringToCSVMolecularMatch.append(stringPriorityTags)
                            .append(";")
                            .append(stringCompositeKeyTags)
                            .append(";")
                            .append(stringSuppressTags)
                            .append(";")
                            .append(stringFilterTypeTags)
                            .append(";")
                            .append(stringTermTags)
                            .append(";")
                            .append(stringPrimaryTags)
                            .append(";")
                            .append(stringFacetTags)
                            .append(";")
                            .append(stringValidTags)
                            .append(";")
                            .append(stringCustomTags)
                            .append(";");

                } else if (keysOfMolecularMatch.get(x).equals("classifications")) {
                    JsonArray molecluarMatchArray =
                            object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchArray.size(); i++) {
                        JsonObject objectclassifications = (JsonObject) molecluarMatchArray.get(i);
                        List<String> keysclassifications = new ArrayList<>(objectclassifications.keySet());
                        for (int t = 0; t < keysclassifications.size(); t++) {
                            stringToCSVMolecularMatch.append(objectclassifications.get(keysclassifications.get(t))).append(";");
                        }
                    }
                } else if (keysOfMolecularMatch.get(x).equals("therapeuticContext")) {
                    JsonArray molecluarMatchArray =
                            object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)).getAsJsonArray();
                    for (int i = 0; i < molecluarMatchArray.size(); i++) {
                        JsonObject objecttherapeuticContext = (JsonObject) molecluarMatchArray.get(i);
                        List<String> keystherapeuticContext = new ArrayList<>(objecttherapeuticContext.keySet());
                        for (int t = 0; t < keystherapeuticContext.size(); t++) {
                            stringToCSVMolecularMatch.append(objecttherapeuticContext.get(keystherapeuticContext.get(t))).append(";");
                        }
                    }
                } else {
                    stringToCSVMolecularMatch.append(object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x))).append(";");
                }
            }
        } else {

        }

        return stringToCSVMolecularMatch;
    }

}
