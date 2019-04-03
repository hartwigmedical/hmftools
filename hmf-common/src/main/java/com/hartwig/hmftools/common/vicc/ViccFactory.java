package com.hartwig.hmftools.common.vicc;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonPrimitive;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public abstract class ViccFactory {
    private static final Logger LOGGER = LogManager.getLogger(ViccFactory.class);

    private ViccFactory() {
    }

    private static StringBuilder readObjectBRCA(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //BRCA object
        StringBuilder stringToCSVBRCA = new StringBuilder();
        if (object.getAsJsonObject("brca") != null) {
            for (int i = 0; i < object.getAsJsonObject("brca").keySet().size(); i++) {
                List<String> keysOfBRCAObject = new ArrayList<>(object.getAsJsonObject("brca").keySet());
                stringToCSVBRCA.append(object.getAsJsonObject("brca").get(keysOfBRCAObject.get(i))).append(";"); // brca data
            }
        } else {
        }

        return stringToCSVBRCA;
    }

    private static StringBuilder readObjectCGI(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //CGI object
        StringBuilder stringToCSVCGI = new StringBuilder();
        StringBuilder stringInfo = new StringBuilder();
        if (object.getAsJsonObject("cgi") != null) {
            List<String> keysOfCGI = new ArrayList<>(object.getAsJsonObject("cgi").keySet());
            for (int i = 0; i < keysOfCGI.size(); i++) {
                if (keysOfCGI.get(i).equals("info")) {
                    JsonElement elementInfo = object.getAsJsonObject("cgi").get(keysOfCGI.get(i));
                    for (int z = 0; z < elementInfo.getAsJsonArray().size(); z++) {
                        String info = elementInfo.getAsJsonArray().get(z).toString().replaceAll(";", ":");
                        stringInfo.append(info).append(",");
                    }
                    stringToCSVCGI.append(stringInfo).append(";");
                } else {
                    stringToCSVCGI.append(object.getAsJsonObject("cgi").get(keysOfCGI.get(i))).append(";");
                }
            }
        } else {
            stringToCSVCGI.append(";;;;;;;;;;;;;;;;;;;;;; ");

        }
        return stringToCSVCGI;
    }

    private static void commandCIVIC(@NotNull JsonObject civicObject, @NotNull List<String> keysOfLifeCycleActions, @NotNull String i,
            @NotNull StringBuilder stringToCSVCIVIC) {
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

    private static StringBuilder readObjectCIVIC(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
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
            LOGGER.info(keysOfCivic);
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
                            commandCIVIC(civicObject, keysOfLifeCycleActions, "last_commented_on", stringToCSVCIVIC);
                        }
                        if (keysOfLifeCycleActions.get(i).equals("last_modified")) {
                            commandCIVIC(civicObject, keysOfLifeCycleActions, "last_modified", stringToCSVCIVIC);
                        }
                        if (keysOfLifeCycleActions.get(i).equals("last_reviewed")) {
                            commandCIVIC(civicObject, keysOfLifeCycleActions, "last_reviewed", stringToCSVCIVIC);
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
        headerCSV.append("variant_groups; entrez_name; display_name;description;url;so_id;id;name; civic_actionability_score; "
                + "clinvar_entries; last_commented_on.timestamp;last_commented_on.username;last_commented_on.area_of_expertise;"
                + "last_commented_on.url;last_commented_on.id;last_commented_on.x32;last_commented_on.x256;last_commented_on.x14;"
                + "last_commented_on.x64;last_commented_on.x128; last_commented_on.description; last_commented_on.name;  "
                + "last_commented_on.twitter_handle;  last_commented_on.name;   last_commented_on.bio;  "
                + "last_commented_on.url; last_commented_on.created_at;  last_commented_on.x32;  last_commented_on.x14; "
                + "last_commented_on.x64; last_commented_on.x128; last_commented_on.accepted_license;last_commented_on.affiliation;"
                + "last_commented_on.avatar_url;last_commented_on.role;last_commented_on.facebook_profile;"
                + "last_commented_on.linkedin_profile; last_commented_on.orcid; last_commented_on.display_name; "
                + "last_commented_on.last_seen_at; last_commented_on.featured_expert; last_commented_on.id;"
                + " last_commented_on.signup_complete; last_modified.timestamp;last_modified.username;last_modified.area_of_expertise;"
                + "last_modified.url;last_modified.id;last_modified.x32;last_modified.x256;last_modified.x14;last_modified.x64;"
                + "last_modified.x128; last_modified.description; last_modified.name;  last_modified.twitter_handle;  last_modified.name;  "
                + " last_modified.bio;  last_modified.url; last_modified.created_at;  last_modified.x32;  last_modified.x14; "
                + "last_modified.x64; last_modified.x128; last_modified.accepted_license;last_modified.affiliation;last_modified.avatar_url;"
                + "last_modified.role;last_modified.facebook_profile;last_modified.linkedin_profile; last_modified.orcid; "
                + "last_modified.display_name; last_modified.last_seen_at; last_modified.featured_expert; last_modified.id; "
                + "last_modified.signup_complete;last_reviewed.timestamp;last_reviewed.username;last_reviewed.area_of_expertise;"
                + "last_reviewed.url;last_reviewed.id;last_reviewed.x32;last_reviewed.x256;last_reviewed.x14;last_reviewed.x64;"
                + "last_reviewed.x128; last_reviewed.description; last_reviewed.name;  last_reviewed.twitter_handle;  last_reviewed.name; "
                + "  last_reviewed.bio;  last_reviewed.url; last_reviewed.created_at;  last_reviewed.x32;  last_reviewed.x14; "
                + "last_reviewed.x64; last_reviewed.x128; last_reviewed.accepted_license;last_reviewed.affiliation;last_reviewed.avatar_url;"
                + "last_reviewed.role;last_reviewed.facebook_profile;last_reviewed.linkedin_profile; last_reviewed.orcid; "
                + "last_reviewed.display_name; last_reviewed.last_seen_at; last_reviewed.featured_expert; last_reviewed.id; "
                + "last_reviewed.signup_complete; variant_aliases; allele_registry_id; provisional_values; gene_id; name; status; "
                + "rating; drug_interaction_type; description; open_change_count; evidence_type; pubchem_id; id; name; variant_origin;"
                + " doid;url;display_name;id;name;status;open_access;name;journal;citation;pmc_id;full_journal_title;source_url;"
                + "clinical_trials;pubmed_id;is_review;year;day;month;id;evidence_direction;variant_id;clinical_significance;evidence_level;"
                + "type;id;name; sources; entrez_id; assertions; hgvs_expressions; errors; chromosome2;reference_bases;start2;variant_bases;"
                + "stop;stop2;representative_transcript2;start;representative_transcript;ensembl_version;chromosome;reference_build; type;"
                + " id; description");
        return stringToCSVCIVIC;
    }

    private static StringBuilder readObjectJax(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //Jax object
        StringBuilder stringToCSVJax = new StringBuilder();
        StringBuilder stringUrl = new StringBuilder();
        StringBuilder stringId = new StringBuilder();
        StringBuilder stringPubmedId = new StringBuilder();
        StringBuilder stringTitle = new StringBuilder();
        if (object.getAsJsonObject("jax") != null) {
            List<String> keysOfJax = new ArrayList<>(object.getAsJsonObject("jax").keySet());
            for (int j = 0; j < keysOfJax.size(); j++) {
                if (keysOfJax.get(j).equals("molecularProfile")) {
                    JsonObject jaxObject = object.getAsJsonObject("jax").get(keysOfJax.get(j)).getAsJsonObject();
                    List<String> keysOfMolecularProfile = new ArrayList<>(jaxObject.keySet());
                    for (int x = 0; x < jaxObject.keySet().size(); x++) {
                        stringToCSVJax.append(jaxObject.get(keysOfMolecularProfile.get(x))).append(";");
                    }
                } else if (keysOfJax.get(j).equals("therapy")) {
                    JsonObject jaxObject = object.getAsJsonObject("jax").get(keysOfJax.get(j)).getAsJsonObject();
                    List<String> keysOfTherapy = new ArrayList<>(jaxObject.keySet());
                    for (int x = 0; x < jaxObject.keySet().size(); x++) {
                        stringToCSVJax.append(jaxObject.get(keysOfTherapy.get(x))).append(";");
                    }
                } else if (keysOfJax.get(j).equals("indication")) {
                    JsonObject jaxObject = object.getAsJsonObject("jax").get(keysOfJax.get(j)).getAsJsonObject();
                    List<String> keysOfIndication = new ArrayList<>(jaxObject.keySet());
                    for (int x = 0; x < jaxObject.keySet().size(); x++) {
                        stringToCSVJax.append(jaxObject.get(keysOfIndication.get(x))).append(";");
                    }
                } else if (keysOfJax.get(j).equals("references")) {
                    JsonArray jaxArray = object.getAsJsonObject("jax").get(keysOfJax.get(j)).getAsJsonArray();
                    for (int x = 0; x < jaxArray.size(); x++) {
                        JsonObject objectRefereces = (JsonObject) jaxArray.get(x);
                        for (int v = 0; v < objectRefereces.keySet().size(); v++) {
                            List<String> keysRefereces = new ArrayList<>(objectRefereces.keySet());
                            if (keysRefereces.get(v).equals("url")) {
                                JsonElement url = objectRefereces.get(keysRefereces.get(v));
                                stringUrl.append(url).append(",");
                            } else if (keysRefereces.get(v).equals("id")) {
                                JsonElement url = objectRefereces.get(keysRefereces.get(v));
                                stringId.append(url).append(",");
                            } else if (keysRefereces.get(v).equals("pubMedId")) {
                                JsonElement url = objectRefereces.get(keysRefereces.get(v));
                                stringPubmedId.append(url).append(",");
                            } else if (keysRefereces.get(v).equals("title")) {
                                JsonElement url = objectRefereces.get(keysRefereces.get(v));
                                String urlAdapted = url.toString().replace(";", ":");
                                stringTitle.append(urlAdapted).append(",");
                            }
                        }
                    }
                    stringToCSVJax.append(stringUrl)
                            .append(";")
                            .append(stringId)
                            .append(";")
                            .append(stringPubmedId)
                            .append(";")
                            .append(stringTitle)
                            .append(";");
                } else {
                    stringToCSVJax.append(object.getAsJsonObject("jax").get(keysOfJax.get(j))).append(";");
                }
            }

        } else {
            stringToCSVJax.append(";;;;;;;;;;;;;;;;");
        }
        return stringToCSVJax;
    }

    private static StringBuilder readObjectJaxTrials(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
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
            stringToCSVJaxTrials.append(";;;;;;;;;;;;;;;");
        }
        return stringToCSVJaxTrials;
    }

    private static StringBuilder readObjectMolecularMatch(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
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
            LOGGER.info(keysOfMolecularMatch);
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
                                for (int p = 0; p < arrayTranscriptConsequence.size(); p++) {
                                    JsonObject objectTranscriptConsequence = (JsonObject) arrayTranscriptConsequence.get(u);
                                    List<String> keysTranscriptConsequence = new ArrayList<>(objectTranscriptConsequence.keySet());
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
                            } else if (keysMutations.get(u).equals("GRCh37_location")) {
                                JsonArray arrayLocationRefGenome = objectMutations.get("GRCh37_location").getAsJsonArray();
                                for (int r = 0; r < arrayLocationRefGenome.size(); r++) {
                                    JsonObject objectLocationRefGenome = (JsonObject) arrayLocationRefGenome.get(i);
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
                    stringToCSVMolecularMatch.append(object.getAsJsonObject("molecularmatch").get(keysOfMolecularMatch.get(x)));
                }
            }

        } else {

        }
        return stringToCSVMolecularMatch;
    }

    private static StringBuilder readObjectMolecularMatchTrials(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
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
            LOGGER.info(object.getAsJsonObject("molecularmatch_trials"));
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
                                stringContact.append(contact).append(",");
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
                                            JsonElement type = objectLocations.get(keysLocations.get(p));
                                            stringType.append(type).append(",");
                                        } else if (keys.get(y).equals("coordinates")) {
                                            JsonElement coordinates = objectLocations.get(keysLocations.get(p));
                                            stringCoordinates.append(coordinates).append(",");
                                        }
                                    }
                                }
                                stringLocation.append(stringType).append(",");
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
                                        JsonElement lat = objectLocations.get(keysLocations.get(p));
                                        stringLat.append(lat).append(",");
                                    } else if (keys.get(v).equals("lon")) {
                                        JsonElement lon = objectLocations.get(keysLocations.get(p));
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
                    stringToCSVMolecularMatchTrials.append(stringInterventionName)
                            .append(";")
                            .append(stringStatus)
                            .append(";")
                            .append(stringCity)
                            .append(";")
                            .append(stringValid)
                            .append(";")
                            .append(stringZip)
                            .append(";")
                            .append(stringCreated)
                            .append(";")
                            .append(stringId)
                            .append(";")
                            .append(stringLastUpdated)
                            .append(";")
                            .append(stringContact)
                            .append(";")
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
                    JsonObject objectOverallContact =
                            object.getAsJsonObject("molecularmatch_trials").get(keysOfMolecularMatchTrials.get(x)).getAsJsonObject();
                    List<String> keysOfOverallContact = new ArrayList<>(objectOverallContact.keySet());
                    for (int u = 0; u < keysOfOverallContact.size(); u++) {
                        stringToCSVMolecularMatchTrials.append(objectOverallContact.get(keysOfOverallContact.get(u)));
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
                    stringToCSVMolecularMatchTrials.append(stringInterventionName)
                            .append(";")
                            .append(stringFacet)
                            .append(";")
                            .append(stringCompositeKey)
                            .append(";")
                            .append(stringSuppress)
                            .append(";")
                            .append(stringGenerateBy)
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
                            .append(";");
                } else {
                    //                    JsonObject molecularMatchTrialsObject =
                    //                            object.getAsJsonObject("molecularmatch_trials").get(keysOfMolecularMatchTrials.get(x)).getAsJsonObject();
                    //                    stringToCSVMolecularMatchTrials.append(molecularMatchTrialsObject.get(keysOfMolecularMatchTrials.get(x)));
                }
            }

        } else {
            headerCSV.append("");
        }
        return stringToCSVMolecularMatchTrials;
    }

    private static StringBuilder readObjectOncokb(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
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
            stringToCSVOncoKb.append(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;");
        }
        return stringToCSVOncoKb;
    }

    private static StringBuilder readObjectPmkb(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //PMKB object
        StringBuilder stringToCSVPmkb = new StringBuilder();
        List<String> keysOfPmkb;
        String header = "";
        int indexValue;
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
                    headerCSV.append(String.join(";", pmkbObject.keySet())).append(";");
                } else if (keysOfPmkb.get(j).equals("tissues")) {
                    JsonArray arrayTissue = object.getAsJsonObject("pmkb").get(keysOfPmkb.get(j)).getAsJsonArray();
                    for (int x = 0; x < arrayTissue.size(); x++) {
                        JsonObject objectTissue = (JsonObject) arrayTissue.get(x);
                        Set<String> set = objectTissue.keySet();
                        header = String.join(";", set);
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
                    headerCSV.append(String.join(";", header)).append(";");
                    stringToCSVPmkb.append(stringId).append(";").append(stringName).append(";");
                } else if (keysOfPmkb.get(j).equals("variant")) {
                    JsonObject pmkbObject = object.getAsJsonObject("pmkb").get(keysOfPmkb.get(j)).getAsJsonObject();
                    List<String> keysOfVariant = new ArrayList<>(pmkbObject.keySet());
                    for (int x = 0; x < pmkbObject.keySet().size(); x++) {
                        if (keysOfVariant.get(x).equals("gene")) {
                            JsonElement elementGene = object.getAsJsonObject("pmkb").get("variant");
                            List<String> keysGene = new ArrayList<>(elementGene.getAsJsonObject().get("gene").getAsJsonObject().keySet());

                            indexValue = keysOfVariant.indexOf("gene");
                            keysOfVariant.remove(indexValue);
                            keysOfVariant.add(indexValue, String.join(";", keysGene));

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

    private static StringBuilder readObjectSage(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
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

    private static StringBuilder readObjectSource(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //Source object
        StringBuilder stringToCSVSource = new StringBuilder();
        stringToCSVSource.append(object.getAsJsonPrimitive("source")).append(";"); // source data
        return stringToCSVSource;
    }

    private static StringBuilder readObjectGenes(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //Genes object
        StringBuilder stringToCSVGenes = new StringBuilder();
        JsonArray arrayGenes = object.getAsJsonArray("genes");
        String genes = arrayGenes.toString();
        genes = genes.substring(1, genes.length() - 1);
        stringToCSVGenes.append(genes).append(";"); // genes data
        return stringToCSVGenes;
    }

    private static StringBuilder readObjectTags(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //Tags object
        StringBuilder stringToCSVTags = new StringBuilder();
        JsonArray arrayTags = object.getAsJsonArray("tags");
        String tags = arrayTags.toString();
        tags = tags.substring(1, tags.length() - 1);
        stringToCSVTags.append(tags).append(";"); // tags data
        return stringToCSVTags;
    }

    private static StringBuilder readObjectDevTags(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //dev_tags
        StringBuilder stringToCSVDevTags = new StringBuilder();
        JsonArray arrayDevTags = object.getAsJsonArray("dev_tags");
        String devTags = arrayDevTags.toString();
        devTags = devTags.substring(1, devTags.length() - 1);
        stringToCSVDevTags.append(devTags).append(";"); // dev tags data
        return stringToCSVDevTags;
    }

    private static StringBuilder readObjectGeneIdentifiers(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //gene_identifiers object
        StringBuilder stringToCSVGeneIdentifiers = new StringBuilder();
        StringBuilder stringSymbol = new StringBuilder();
        StringBuilder stringEntrezId = new StringBuilder();
        StringBuilder stringEnsembleGeneId = new StringBuilder();
        JsonArray arrayGeneIdentifiers = object.getAsJsonArray("gene_identifiers");
        if (arrayGeneIdentifiers.size() == 0) {
            stringToCSVGeneIdentifiers.append(";;;");
        } else {
            for (int j = 0; j < arrayGeneIdentifiers.size(); j++) {
                JsonObject objectGeneIdentiefiers = (JsonObject) arrayGeneIdentifiers.get(j);
                for (int i = 0; i < objectGeneIdentiefiers.keySet().size(); i++) {
                    List<String> keysOfGeneIdentifiersObject = new ArrayList<>(objectGeneIdentiefiers.keySet());
                    if (keysOfGeneIdentifiersObject.get(i).equals("symbol")) {
                        JsonElement symbol = objectGeneIdentiefiers.get(keysOfGeneIdentifiersObject.get(i));
                        stringSymbol.append(symbol).append(",");
                    } else if (keysOfGeneIdentifiersObject.get(i).equals("entrez_id")) {
                        JsonElement entrezId = objectGeneIdentiefiers.get(keysOfGeneIdentifiersObject.get(i));
                        stringEntrezId.append(entrezId).append(",");
                    } else if (keysOfGeneIdentifiersObject.get(i).equals("ensembl_gene_id")) {
                        JsonElement ensemblGeneID = objectGeneIdentiefiers.get(keysOfGeneIdentifiersObject.get(i));
                        stringEnsembleGeneId.append(ensemblGeneID).append(",");
                    }
                }
            }
            stringToCSVGeneIdentifiers.append(stringSymbol)
                    .append(";")
                    .append(stringEntrezId)
                    .append(";")
                    .append(stringEnsembleGeneId)
                    .append(";");
        }
        return stringToCSVGeneIdentifiers;
    }

    private static StringBuilder readObjectAssociation(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //association object
        List<String> keysOfAssocationObject = Lists.newArrayList();

        StringBuilder stringToCSVAssociation = new StringBuilder();
        for (int i = 0; i < object.getAsJsonObject("association").keySet().size(); i++) {
            keysOfAssocationObject = new ArrayList<>(object.getAsJsonObject("association").keySet());
            LOGGER.info(keysOfAssocationObject);

            if (keysOfAssocationObject.get(i).equals("description")) {
                stringToCSVAssociation.append(object.getAsJsonObject("association").get(keysOfAssocationObject.get(i)))
                        .append(";"); // association data
            } else if (keysOfAssocationObject.get(i).equals("evidence")) {
                JsonElement elementEvidence = object.getAsJsonObject("association").get("evidence");
                JsonArray arrayEvidence = elementEvidence.getAsJsonArray();
                JsonObject objectEvidence = (JsonObject) arrayEvidence.iterator().next();
                LOGGER.info(objectEvidence.keySet());

                for (int a = 0; a < objectEvidence.keySet().size(); a++) {
                    List<String> keysOfEvidenceObject = new ArrayList<>(objectEvidence.keySet());
                    if (keysOfEvidenceObject.get(a).equals("evidenceType")) {
                        for (int b = 0; b < objectEvidence.get("evidenceType").getAsJsonObject().keySet().size(); b++) {
                            List<String> keysOfEvidenceTypeObject =
                                    new ArrayList<>(objectEvidence.get("evidenceType").getAsJsonObject().keySet());
                            stringToCSVAssociation.append(objectEvidence.get("evidenceType")
                                    .getAsJsonObject()
                                    .get(keysOfEvidenceTypeObject.get(b))).append(";"); // association data
                        }
                    } else {
                        stringToCSVAssociation.append(objectEvidence.get(keysOfEvidenceObject.get(a))).append(";"); // association data
                    }
                }
            } else if (keysOfAssocationObject.get(i).equals("environmentalContexts")) {
                stringToCSVAssociation.append(object.getAsJsonObject("association").get("environmentalContexts"))
                        .append(";"); // association data
            } else if (keysOfAssocationObject.get(i).equals("evidence_label")) {
                stringToCSVAssociation.append(object.getAsJsonObject("association").get(keysOfAssocationObject.get(i)))
                        .append(";"); // association data
            } else if (keysOfAssocationObject.get(i).equals("phenotype")) {
                JsonElement elementPhenotype = object.getAsJsonObject("association").get("phenotype");
                for (int a = 0; a < elementPhenotype.getAsJsonObject().keySet().size(); a++) {
                    List<String> keysOfPhenotypeObject = new ArrayList<>(elementPhenotype.getAsJsonObject().keySet());
                    if (keysOfPhenotypeObject.get(a).equals("type")) {

                        List<String> keysOfPhenotypeTypeObject =
                                new ArrayList<>(elementPhenotype.getAsJsonObject().get("type").getAsJsonObject().keySet());

                        for (int c = 0; c < keysOfPhenotypeObject.size(); c++) {
                            stringToCSVAssociation.append(elementPhenotype.getAsJsonObject() // association data
                                    .get("type").getAsJsonObject().get(keysOfPhenotypeTypeObject.get(c))).append(";");
                        }
                    } else {
                        stringToCSVAssociation.append(elementPhenotype.getAsJsonObject().get(keysOfPhenotypeObject.get(a)))
                                .append(";"); // association data
                    }
                }
            } else if (keysOfAssocationObject.get(i).equals("oncogenic")) {
                stringToCSVAssociation.append(object.getAsJsonObject("association").get(keysOfAssocationObject.get(i)))
                        .append(";"); // association data
            }
        }
        headerCSV.append(keysOfAssocationObject).append(";"); // header features
        return stringToCSVAssociation;
    }

    private static StringBuilder readObjectFeaturesNames(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //feature_names object
        StringBuilder stringToCSVFeaturesNames = new StringBuilder();
        stringToCSVFeaturesNames.append(object.getAsJsonPrimitive("feature_names")).append(";"); // features names data
        headerCSV.append("feature_names").append(";"); // header features names
        return stringToCSVFeaturesNames;
    }

    private static StringBuilder readObjectFeatures(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //features
        StringBuilder stringToCSVFeatures = new StringBuilder();
        List<String> b = Lists.newArrayList();
        List<String> keysOfSequenceOntologyObjectMerged;
        JsonArray arrayFeatures = object.getAsJsonArray("features");
        if (arrayFeatures.size() == 0) {
            stringToCSVFeatures.append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";")
                    .append(Strings.EMPTY)
                    .append(";");

            headerCSV.append("provenance_rule")
                    .append(";")
                    .append("end")
                    .append(";")
                    .append("description")
                    .append(";")
                    .append("links")
                    .append(";")
                    .append("hierarchy")
                    .append(";")
                    .append("soid")
                    .append(";")
                    .append("parent_soid")
                    .append(";")
                    .append("name")
                    .append(";")
                    .append("parent_name")
                    .append(";")
                    .append("provenance")
                    .append(";")
                    .append("start")
                    .append(";")
                    .append("synonyms")
                    .append(";")
                    .append("biomarker_type")
                    .append(";")
                    .append("referenceName")
                    .append(";")
                    .append("alt")
                    .append(";")
                    .append("ref")
                    .append(";")
                    .append("chromosome")
                    .append(";")
                    .append("name")
                    .append(";");
        } else {
            JsonObject objectFeatures = (JsonObject) arrayFeatures.iterator().next();
            int indexValue;
            List<String> keysOfFeaturesObject = new ArrayList<>(objectFeatures.keySet());

            for (int i = 0; i < objectFeatures.keySet().size(); i++) {

                if (keysOfFeaturesObject.get(i).equals("sequence_ontology")) {
                    List<String> keysOfSequenceOntologyObject =
                            new ArrayList<>(objectFeatures.getAsJsonObject("sequence_ontology").keySet());
                    indexValue = keysOfFeaturesObject.indexOf("sequence_ontology");
                    keysOfFeaturesObject.set(indexValue, String.join(";", keysOfSequenceOntologyObject));

                    for (int e = 0; e < objectFeatures.getAsJsonObject("sequence_ontology").keySet().size(); e++) {
                        b.add(objectFeatures.get("sequence_ontology")
                                .getAsJsonObject()
                                .get(keysOfSequenceOntologyObject.get(e))
                                .toString());

                    }
                } else {
                    b.add(objectFeatures.get(keysOfFeaturesObject.get(i)).toString());
                }

            }

            keysOfSequenceOntologyObjectMerged = keysOfFeaturesObject;
            if (!keysOfFeaturesObject.get(1).equals("entrez_id")) {
                keysOfSequenceOntologyObjectMerged.add(1, "entrez_id");
                b.add(1, Strings.EMPTY);
            }
            if (!keysOfFeaturesObject.get(3).equals("name")) {
                keysOfSequenceOntologyObjectMerged.add(3, "name");
                b.add(3, Strings.EMPTY);
            }
            stringToCSVFeatures.append(b.toString().replace(",", ";"));

            String header = String.join(";", keysOfSequenceOntologyObjectMerged);
            headerCSV.append(header).append(";"); // header features
        }

        return stringToCSVFeatures;
    }

    public static void extractAllFile(@NotNull String allJsonPath) throws IOException {
        final String csvFileName = "/Users/liekeschoenmaker/hmf/tmp/all.csv";
        PrintWriter writer = new PrintWriter(new File(csvFileName));
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(allJsonPath));
        reader.setLenient(true);
        int index = 1;
        StringBuilder headerCSV = new StringBuilder();

        String headerIndex = "index;";
        String headerSource = "source;";
        String headerGenes = "genes;";
        String headerTags = "tags;";
        String headerDevTags = "dev_tags;";
        String headerGeneIdentifiers = "gene_identifiers.Symbol;gene_identifiers.entrez_id;gene_identifiers.ensembl_gene_id;";
        String headerSage = "sage.entrez_id;sage.clinical_manifestation;sage.publication_url;sage.germline_or_somatic;sage.evidence_label;"
                + "sage.drug_labels;sage.response_type;sage.gene;";
        String headerPmkb = "pmkb.tumor.id;pmkb.tumor.name;pmkb.tissues.id;pmkb.tissues.name;pmkb.variant.amino_acid_change;"
                + "pmkb.variant.germline;pmkb.variant.partner_gene;pmkb.variant.codons;pmkb.variant.description;pmkb.variant.exons;"
                + "pmkb.variant.notes;pmkb.variant.cosmic;pmkb.variant.effect;pmkb.variant.cnv_type;pmkb.variant.id;"
                + "pmkb.variant.cytoband;pmkb.variant.variant_type;pmkb.variant.dna_change;pmkb.variant.coordinates;"
                + "pmkb.variant.chromosome_based_cnv;pmkb.variant.gene.description;pmkb.variant.gene.created_at;"
                + "pmkb.variant.gene.updated_at;pmkb.variant.gene.active_ind;pmkb.variant.gene.external_id;"
                + "pmkb.variant.gene.id;pmkb.variant.gene.name;pmkb.transcript;pmkb.description_type;pmkb.chromosome;pmkb.name;";
        String headerBRCA = "brca.Variant_frequency_LOVD;brca.Allele_frequency_FIN_ExAC;brca.ClinVarAccession_ENIGMA;"
                + "brca.Homozygous_count_AFR_ExAC;brca.BX_ID_ExAC;brca.Variant_in_LOVD;brca.Allele_frequency_AFR_ExAC;brca.DBID_LOVD;"
                + "brca.Chr;brca.BX_ID_ENIGMA;brca.Co_occurrence_LR_exLOVD;brca.Homozygous_count_EAS_ExAC;brca.Submitter_ClinVar;"
                + "brca.Allele_frequency_EAS_ExAC;brca.Hg37_End;brca.Submitters_LOVD;brca.Clinical_classification_BIC;"
                + "brca.Homozygous_count_NFE_ExAC;brca.Allele_count_SAS_ExAC;brca.Method_ClinVar;brca.Allele_count_NFE_ExAC;"
                + "brca.Pathogenicity_all;brca.Germline_or_Somatic_BIC;brca.Homozygous_count_SAS_ExAC;brca.BIC_Nomenclature;"
                + "brca.Assertion_method_ENIGMA;brca.Literature_source_exLOVD;brca.Change_Type_id;brca.Collection_method_ENIGMA;"
                + "brca.Sum_family_LR_exLOVD;brca.HGVS_cDNA_LOVD”;brca.Homozygous_count_FIN_ExAC;brca.EAS_Allele_frequency_1000_Genomes;"
                + "brca.Ethnicity_BIC;brca.Individuals_LOVD;brca.Variant_in_ExAC;brca.URL_ENIGMA;brca.Allele_Origin_ClinVar;"
                + "brca.Allele_frequency_AMR_ExAC;brca.Variant_in_1000_Genomes;brca.AFR_Allele_frequency_1000_Genomes;"
                + "brca.BX_ID_exLOVD;brca.Source;brca.Condition_ID_value_ENIGMA;brca.HGVS_Protein;brca.Ref;brca.Allele_number_AFR_ExAC;"
                + "brca.Allele_count_AFR_ExAC;brca.BX_ID_LOVD;brca.Synonyms;brca.Gene_Symbol;brca.Comment_on_clinical_significance_ENIGMA;"
                + "brca.Missense_analysis_prior_probability_exLOVD;brca.Allele_number_FIN_ExAC;brca.Posterior_probability_exLOVD;"
                + "brca.Polyphen_Score;brca.Reference_Sequence;brca.Allele_count_EAS_ExAC;brca.Hg38_End;brca.HGVS_cDNA;"
                + "brca.Functional_analysis_technique_LOVD;brca.SAS_Allele_frequency_1000_Genomes;brca.RNA_LOVD;"
                + "brca.Combined_prior_probablility_exLOVD;brca.BX_ID_ClinVar;brca.IARC_class_exLOVD;brca.BX_ID_BIC;brca.Sift_Prediction;"
                + "brca.Allele_number_NFE_ExAC;brca.Allele_origin_ENIGMA”;brca.Allele_number_OTH_ExAC;brca.Hg36_End;"
                + "brca.Allele_frequency_SAS_ExAC;brca.Date_Last_Updated_ClinVar;brca.Allele_number_EAS_ExAC;"
                + "brca.Allele_frequency_OTH_ExAC;brca.Source_URL;brca.SCV_ClinVar;brca.Pathogenicity_expert;"
                + "brca.Allele_frequency_1000_Genomes;brca.Functional_analysis_result_LOVD;brca.AMR_Allele_frequency_1000_Genomes;"
                + "brca.Variant_in_ESP;brca.Variant_in_BIC;brca.Clinical_significance_ENIGMA;brca.Max_Allele_Frequency;"
                + "brca.Allele_count_AMR_ExAC;brca.Variant_in_ENIGMA;brca.BX_ID_ESP;brca.Patient_nationality_BIC”;brca.BX_ID_1000_Genomes;"
                + "brca.Genomic_Coordinate_hg37;brca.Genomic_Coordinate_hg36;brca.EUR_Allele_frequency_1000_Genomes;"
                + "brca.Number_of_family_member_carrying_mutation_BIC;brca.Segregation_LR_exLOVD”;brca.Allele_Frequency;"
                + "brca.Minor_allele_frequency_percent_ESP;brca.Allele_frequency_ExAC”;brca.Mutation_type_BIC;"
                + "brca.Assertion_method_citation_ENIGMA;brca.Condition_ID_type_ENIGMA;brca.Allele_count_OTH_ExAC;brca.HGVS_protein_LOVD;"
                + "brca.Variant_in_ClinVar;brca.Clinical_importance_BIC;brca.Discordant;brca.Allele_count_FIN_ExAC;"
                + "brca.Condition_category_ENIGMA;brca.Allele_Frequency_ESP;brca.Homozygous_count_OTH_ExAC;brca.Genetic_origin_LOVD;"
                + "brca.id;brca.Homozygous_count_AMR_ExAC;brca.Clinical_Significance_ClinVar;brca.AA_Allele_Frequency_ESP;"
                + "brca.Protein_Change;brca.Variant_in_exLOVD;brca.EA_Allele_Frequency_ESP;brca.HGVS_RNA;"
                + "brca.Clinical_significance_citations_ENIGMA;brca.Variant_effect_LOVD;brca.Polyphen_Prediction;brca.Data_Release_id;"
                + "brca.Hg37_Start;brca.Hg36_Start;brca.Sift_Score;brca.Genomic_Coordinate_hg38;brca.Alt;brca.Literature_citation_BIC;"
                + "brca.Variant_haplotype_LOVD;brca.Allele_frequency_NFE_ExAC;brca.Hg38_Start;brca.Pos;brca.Date_last_evaluated_ENIGMA;"
                + "brca.Allele_number_SAS_ExAC;brca.Allele_number_AMR_ExAC;";
        String headerCGI = "cgi.Targeting;cgi.Source”;cgi.cDNA;cgi.Primary Tumor type;cgi.individual_mutation;cgi.Drug full name;"
                + "cgi.Curator;cgi.Drug family;cgi.Alteration;cgi.Drug;cgi.Biomarker;cgi.gDNA;cgi.Drug status;cgi.Gene;cgi.transcript;"
                + "cgi.strand;cgi.info;cgi.Assay type;cgi.Alteration type;cgi.region;cgi.Evidence level;cgi.Association;"
                + "cgi.Metastatic Tumor Type;";
        String headerOncokb = "oncokb.biological.mutationEffectPmids;oncokb.biological.Isoform;oncokb.clinical.RefSeq;"
                + "oncokb.clinical.level;oncokb.clinical.Isoform;oncokb.biological.variant.variantResidues;"
                + "oncokb.biological.variant.proteinStart;oncokb.biological.variant.name;oncokb.biological.variant.proteinEnd;"
                + "oncokb.biological.variant.refResidues;oncokb.biological.variant.alteration;oncokb.biological.variant.consequence.term;"
                + "oncokb.biological.variant.consequence.description;oncokb.biological.variant.consequence.isGenerallyTruncating;"
                + "oncokb.biological.variant.gene.oncogene;oncokb.biological.variant.gene.name;oncokb.biological.variant.gene.hugoSymbol;"
                + "oncokb.biological.variant.gene.curatedRefSeq;oncokb.biological.variant.gene.entrezGeneId;"
                + "oncokb.biological.variant.gene.geneAliases;oncokb.biological.variant.gene.tsg;"
                + "oncokb.biological.variant.gene.curatedIsoform;oncokb.biological.Entrez Gene ID;oncokb.biological.oncogenic;"
                + "oncokb.biological.mutationEffect;oncokb.biological.RefSeq;oncokb.biological.gene;"
                + "oncokb.biological.mutationEffectAbstracts;oncokb.clinical.Entrez Gene ID;oncokb.clinical.drugPmids;"
                + "oncokb.clinical.cancerType;oncokb.clinical.drug;oncokb.clinical.drugAbstracts.text;oncokb.clinical.drugAbstracts.link;"
                + "oncokb.clinical.gene;oncokb.clinical.level_label;";
        String headerJax = "jax.responseType;jax.approvalStatus;jax.molecularProfile.profileName;jax.molecularProfile.id;jax.therapy.id;"
                + "jax.therapy.therapyName;jax.evidenceType;jax.indication.source;jax.indication.id;jax.indication.name;"
                + "jax.efficacyEvidence;jax.references.url;jax.references.id;jax.references.pubMedId”;jax.references.title;jax.id;";
        String headerJaxTrials = "jax_trials.indications.source;jax_trials.indications.id;jax_trials.indications.name;jax_trials.title;"
                + "jax_trials.gender;jax_trials.nctId;jax_trials.sponsors;jax_trials.recruitment;jax_trials.variantRequirements;"
                + "jax_trials.updateDate;jax_trials.phase;jax_trials.variantRequirementDetails.molecularProfile.profileName;"
                + "jax_trials.variantRequirementDetails.molecularProfile.id;jax_trials.variantRequirementDetails.requirementType;"
                + "jax_trials.therapies.id;jax_trials.therapies.therapyName;";
        String headerAssociation = "association.drug_labels;association.description;association.publication_url;association.source_link;"
                + "association.variant_name;association.evidence.info.publications;association.evidence.evidenceType.sourceName;"
                + "association.evidence.evidenceType.id;association.evidence.description;association.environmentalContexts.term;"
                + "association.environmentalContexts.description;association.environmentalContexts.taxonomy.kingdom;"
                + "association.environmentalContexts.taxonomy.direct-parent;association.environmentalContexts.taxonomy.class;"
                + "association.environmentalContexts.taxonomy.subclass;association.environmentalContexts.taxonomy.superclass;"
                + "association.environmentalContexts.source;association.environmentalContexts.usan_stem;"
                + "association.environmentalContexts.approved_countries;association.environmentalContexts.id;association.evidence_label;"
                + "association.phenotype.type.source;association.phenotype.type.term;association.phenotype.type.id;"
                + "association.phenotype.description;association.phenotype.family;association.phenotype.id;association.evidence_level;"
                + "association..response_type;";

        headerCSV.append(headerIndex);
        headerCSV.append(headerSource);
        //        headerCSV.append(headerGenes);
        //        headerCSV.append(headerTags);
        //        headerCSV.append(headerDevTags);
        //        headerCSV.append(headerGeneIdentifiers);
        //        headerCSV.append(headerSage);
        //        headerCSV.append(headerPmkb);
        //        headerCSV.append(headerBRCA);
        //        headerCSV.append(headerCGI);
        //        headerCSV.append(headerOncokb);
        //        headerCSV.append(headerJax);
        //        headerCSV.append(headerJaxTrials);
        headerCSV.append(headerAssociation);

        writer.append(headerCSV);
        writer.append("\n");
        while (reader.peek() != JsonToken.END_DOCUMENT && index < 10) {
            LOGGER.info(index);
            JsonObject object = parser.parse(reader).getAsJsonObject();

            StringBuilder stringToCSVAll = new StringBuilder();
            headerCSV.append("index").append(";");
            StringBuilder StringToCSVSource = readObjectSource(object, headerCSV);
            StringBuilder StringToCSVGenes = readObjectGenes(object, headerCSV);
            StringBuilder StringToCSVTags = readObjectTags(object, headerCSV);
            StringBuilder StringToCSVDevTags = readObjectDevTags(object, headerCSV);
            StringBuilder StringToCSVGeneIdentifiers = readObjectGeneIdentifiers(object, headerCSV);
            StringBuilder StringToCSVSage = readObjectSage(object, headerCSV);
            StringBuilder StringToCSVPmkb = readObjectPmkb(object, headerCSV);
            StringBuilder StringToCSVBrca = readObjectBRCA(object, headerCSV);
            StringBuilder StringToCSVCGI = readObjectCGI(object, headerCSV);
            StringBuilder StringToCSVOncokb = readObjectOncokb(object, headerCSV);
            StringBuilder StringToCSVJax = readObjectJax(object, headerCSV);
            StringBuilder StringToCSVJaxTrials = readObjectJaxTrials(object, headerCSV);
          //  StringBuilder StringToCSVAssociation = readObjectAssociation(object, headerCSV);

            stringToCSVAll.append(index).append(";");
            stringToCSVAll.append(StringToCSVSource);
            //            stringToCSVAll.append(StringToCSVGenes);
            //            stringToCSVAll.append(StringToCSVTags);
            //            stringToCSVAll.append(StringToCSVDevTags);
            //            stringToCSVAll.append(StringToCSVGeneIdentifiers);
            //            stringToCSVAll.append(StringToCSVSage);
            //            stringToCSVAll.append(StringToCSVPmkb);
            //            stringToCSVAll.append(StringToCSVBrca);
            //            stringToCSVAll.append(StringToCSVCGI);
            //            stringToCSVAll.append(StringToCSVOncokb);
            //            stringToCSVAll.append(StringToCSVJax);
            //            stringToCSVAll.append(StringToCSVJaxTrials);
          //  stringToCSVAll.append(StringToCSVAssociation);

            writer.append(stringToCSVAll);
            writer.append("\n");
            index++;

        }
        reader.close();
        writer.close();
    }

    public static void extractBRCAFile(@NotNull String brcaJsonPath) throws IOException {
    }

    public static void extractCgiFile(@NotNull String cgiJsonPath) throws IOException {
    }

    public static void extractCivicFile(@NotNull String civicJsonPath) throws IOException {
    }

    public static void extractJaxFile(@NotNull String jaxJsonPath) throws IOException {
    }

    public static void extractJaxTrialsFile(@NotNull String jaxTrialsJsonPath) throws IOException {
    }

    public static void extractMolecularMatchFile(@NotNull String molecularMatchJsonPath) throws IOException {
    }

    public static void extractMolecularMatchTrailsFile(@NotNull String molecularMatchTrialsJsonPath) throws IOException {
    }

    public static void extractOncokbFile(@NotNull String oncokbJsonPath) throws IOException {
    }

    public static void extractPmkbFile(@NotNull String pmkbJsonPath) throws IOException {
    }

    public static void extractSageFile(@NotNull String sageJsonPath) throws IOException {
    }

}
