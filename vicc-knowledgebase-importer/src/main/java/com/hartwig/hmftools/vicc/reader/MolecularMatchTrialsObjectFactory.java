package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.jsonArrayToStringList;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsContact;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsGeo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsIntervation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsLocations;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsOverallContact;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsTags;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsContact;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsGeo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsIntervation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsLocations;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsOverallContact;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsTags;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class MolecularMatchTrialsObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(MolecularMatchTrialsObjectFactory.class);

    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_ELEMENT_SIZES = Lists.newArrayList(13, 14);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_INTERVATIONS_ELEMENT_SIZES = Lists.newArrayList(1, 2, 3, 4, 5);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_LOCATIONS_ELEMENT_SIZES =
            Lists.newArrayList(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_OVERALL_CONTACT_ELEMENT_SIZES =
            Lists.newArrayList(0, 1, 2, 3, 4, 5, 7, 8, 9, 10);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_TAGS_ELEMENT_SIZES = Lists.newArrayList(7, 8, 9, 10, 11, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_GEO_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_LOCATION_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_CONTACT_ELEMENT_SIZES = Lists.newArrayList(0, 1, 2, 3);

    private MolecularMatchTrialsObjectFactory() {
    }

    @NotNull
    static MolecularMatchTrials create(@NotNull JsonObject objectMolecularMatchTrials) {
        Set<String> keysMolecularMatchTrials = objectMolecularMatchTrials.keySet();
        if (!EXPECTED_MOLECULARMATCH_TRAILS_ELEMENT_SIZES.contains(keysMolecularMatchTrials.size())) {
            LOGGER.warn("Found {} in molecular match trials rather than the expected {}",
                    keysMolecularMatchTrials.size(),
                    EXPECTED_MOLECULARMATCH_TRAILS_ELEMENT_SIZES);
            LOGGER.warn(keysMolecularMatchTrials);
        }

        return ImmutableMolecularMatchTrials.builder()
                .status(objectMolecularMatchTrials.getAsJsonPrimitive("status").getAsString())
                .startDate(!objectMolecularMatchTrials.has("startDate") || objectMolecularMatchTrials.get("startDate").isJsonNull()
                        ? null
                        : objectMolecularMatchTrials.getAsJsonPrimitive("startDate").getAsString())
                .title(objectMolecularMatchTrials.getAsJsonPrimitive("title").getAsString())
                .molecularAlterations(jsonArrayToStringList(objectMolecularMatchTrials.getAsJsonArray("molecularAlterations")))
                .score(objectMolecularMatchTrials.getAsJsonPrimitive("_score").getAsString())
                .intervation(createMolecularMatchTrialsInterventions(objectMolecularMatchTrials.getAsJsonArray("interventions")))
                .locations(createMolecularMatchTrialsLocations(objectMolecularMatchTrials.getAsJsonArray("locations")))
                .briefTitle(objectMolecularMatchTrials.get("briefTitle").isJsonNull()
                        ? null
                        : objectMolecularMatchTrials.getAsJsonPrimitive("briefTitle").getAsString())
                .overallContact(objectMolecularMatchTrials.get("overallContact").isJsonNull()
                        ? null
                        : createMolecularMatchTrialsOverallContact(objectMolecularMatchTrials.getAsJsonObject("overallContact")))
                .link(objectMolecularMatchTrials.getAsJsonPrimitive("link").getAsString())
                .phase(objectMolecularMatchTrials.getAsJsonPrimitive("phase").getAsString())
                .tags(createMolecularMatchTrialsTags(objectMolecularMatchTrials.getAsJsonArray("tags")))
                .id(objectMolecularMatchTrials.getAsJsonPrimitive("id").getAsString())
                .studyType(objectMolecularMatchTrials.getAsJsonPrimitive("studyType").getAsString())
                .build();
    }

    @NotNull
    private static List<MolecularMatchTrialsLocations> createMolecularMatchTrialsLocations(@NotNull JsonArray arrayLocations) {
        List<MolecularMatchTrialsLocations> locationsList = Lists.newArrayList();
        for (JsonElement location : arrayLocations) {
            Set<String> keysLocation = location.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRAILS_LOCATIONS_ELEMENT_SIZES.contains(keysLocation.size())) {
                LOGGER.warn("Found {} in molecular match trials locations rather than the expected {}",
                        keysLocation.size(),
                        EXPECTED_MOLECULARMATCH_TRAILS_LOCATIONS_ELEMENT_SIZES);
                LOGGER.warn(keysLocation);
            }

            locationsList.add(ImmutableMolecularMatchTrialsLocations.builder()
                    .status(location.getAsJsonObject().getAsJsonPrimitive("status").getAsString())
                    .last_name(!location.getAsJsonObject().has("last_name")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("last_name").getAsString())
                    .email(!location.getAsJsonObject().has("email")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("email").getAsString())
                    .phone(!location.getAsJsonObject().has("phone")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("phone").getAsString())
                    .phone_backup(!location.getAsJsonObject().has("phone_backup")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("phone_backup").getAsString())
                    .email_backup(!location.getAsJsonObject().has("email_backup")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("email_backup").getAsString())
                    .last_name_backup(!location.getAsJsonObject().has("last_name_backup")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("last_name_backup").getAsString())
                    .phone_ext_backup(!location.getAsJsonObject().has("phone_ext_backup")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("phone_ext_backup").getAsString())
                    .phone_ext(!location.getAsJsonObject().has("phone_ext")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("phone_ext").getAsString())
                    .city(!location.getAsJsonObject().has("city")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("city").getAsString())
                    .valid(!location.getAsJsonObject().has("_valid")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("_valid").getAsString())
                    .zip(!location.getAsJsonObject().has("zip") ? null : location.getAsJsonObject().getAsJsonPrimitive("zip").getAsString())
                    .created(!location.getAsJsonObject().has("created")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("created").getAsString())
                    .country(!location.getAsJsonObject().has("country")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("country").getAsString())
                    .number(!location.getAsJsonObject().has("number")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("number").getAsString())
                    .id(!location.getAsJsonObject().has("id") ? null : location.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .lastUpdated(!location.getAsJsonObject().has("lastUpdated")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("lastUpdated").getAsString())
                    .contact(!location.getAsJsonObject().has("contact")
                            ? null
                            : createMolecularMatchTrialsContact(location.getAsJsonObject().getAsJsonObject("contact")))
                    .state(!location.getAsJsonObject().has("state")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("state").getAsString())
                    .street(!location.getAsJsonObject().has("street")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("street").getAsString())
                    .location(!location.getAsJsonObject().has("location") || location.getAsJsonObject().get("location").isJsonNull()
                            ? null
                            : createMolecularMatchTrialsLocation(location.getAsJsonObject().getAsJsonObject("location")))
                    .po_box(!location.getAsJsonObject().has("po_box")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("po_box").getAsString())
                    .failedGeocode(!location.getAsJsonObject().has("failedGeocode")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("failedGeocode").getAsString())
                    .geo(!location.getAsJsonObject().has("geo")
                            ? null
                            : createMolecularMatchTrialsGeo(location.getAsJsonObject().getAsJsonObject("geo")))
                    .validMessage(!location.getAsJsonObject().has("_validMessage")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("_validMessage").getAsString())
                    .name(!location.getAsJsonObject().has("name")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return locationsList;
    }

    @NotNull
    private static List<MolecularMatchTrialsIntervation> createMolecularMatchTrialsInterventions(@NotNull JsonArray interventionsArray) {
        List<MolecularMatchTrialsIntervation> molecularMatchTrialsInterventionList = Lists.newArrayList();
        for (JsonElement intervention : interventionsArray) {
            Set<String> keysIntervention = intervention.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRAILS_INTERVATIONS_ELEMENT_SIZES.contains(keysIntervention.size())) {
                LOGGER.warn("Found {} in molecular match trials intervention rather than the expected {}",
                        keysIntervention.size(),
                        EXPECTED_MOLECULARMATCH_TRAILS_INTERVATIONS_ELEMENT_SIZES);
                LOGGER.warn(keysIntervention);
            }

            molecularMatchTrialsInterventionList.add(ImmutableMolecularMatchTrialsIntervation.builder()
                    .intervention_name(!intervention.getAsJsonObject().has("intervention_name")
                            ? null
                            : intervention.getAsJsonObject().getAsJsonPrimitive("intervention_name").getAsString())
                    .other_name(!intervention.getAsJsonObject().has("other_name")
                            ? null
                            : otherNameMolecularMatchTrials(intervention.getAsJsonObject()))
                    .description(!intervention.getAsJsonObject().has("description")
                            ? null
                            : intervention.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .arm_group_label(!intervention.getAsJsonObject().has("arm_group_label")
                            ? null
                            : armGroupLabelMolecularMatchTrials(intervention.getAsJsonObject()))
                    .intervention_type(!intervention.getAsJsonObject().has("intervention_type")
                            ? null
                            : intervention.getAsJsonObject().getAsJsonPrimitive("intervention_type").getAsString())
                    .build());
        }
        return molecularMatchTrialsInterventionList;
    }

    @NotNull
    private static Iterable<String> otherNameMolecularMatchTrials(@NotNull JsonObject otherNameObject) {
        if (otherNameObject.get("other_name").isJsonPrimitive()) {
            return Collections.singletonList(otherNameObject.getAsJsonPrimitive("other_name").getAsString());
        } else {
            return jsonArrayToStringList(otherNameObject.getAsJsonArray("other_name"));
        }
    }

    @NotNull
    private static Iterable<String> armGroupLabelMolecularMatchTrials(@NotNull JsonObject armGroupLabel) {
        if (armGroupLabel.get("arm_group_label").isJsonArray()) {
            return jsonArrayToStringList(armGroupLabel.getAsJsonArray("arm_group_label"));
        } else if (armGroupLabel.get("arm_group_label").isJsonPrimitive()) {
            return Collections.singletonList(armGroupLabel.getAsJsonPrimitive("arm_group_label").getAsString());
        } else {
            return Lists.newArrayList();
        }
    }

    @NotNull
    private static MolecularMatchTrialsOverallContact createMolecularMatchTrialsOverallContact(@NotNull JsonObject overallContactObject) {
        Set<String> keysOverallContact = overallContactObject.getAsJsonObject().keySet();
        if (!EXPECTED_MOLECULARMATCH_TRAILS_OVERALL_CONTACT_ELEMENT_SIZES.contains(keysOverallContact.size())) {
            LOGGER.warn("Found {} in molecular match trials overall contact rather than the expected {}",
                    keysOverallContact.size(),
                    EXPECTED_MOLECULARMATCH_TRAILS_OVERALL_CONTACT_ELEMENT_SIZES);
            LOGGER.warn(keysOverallContact);
        }

        return ImmutableMolecularMatchTrialsOverallContact.builder()
                .phone(!overallContactObject.has("phone") || overallContactObject.get("phone").isJsonNull()
                        ? null
                        : overallContactObject.getAsJsonPrimitive("phone").getAsString())
                .last_name(!overallContactObject.has("last_name")
                        ? null
                        : overallContactObject.getAsJsonPrimitive("last_name").getAsString())
                .email(!overallContactObject.has("email") ? null : overallContactObject.getAsJsonPrimitive("email").getAsString())
                .affiliation(!overallContactObject.has("affiliation") || overallContactObject.get("affiliation").isJsonNull()
                        ? null
                        : overallContactObject.getAsJsonPrimitive("affiliation").getAsString())
                .phone_ext(!overallContactObject.has("phone_ext")
                        ? null
                        : overallContactObject.getAsJsonPrimitive("phone_ext").getAsString())
                .country(!overallContactObject.has("country") ? null : overallContactObject.getAsJsonPrimitive("country").getAsString())
                .city(!overallContactObject.has("city") ? null : overallContactObject.getAsJsonPrimitive("city").getAsString())
                .name(!overallContactObject.has("name") ? null : overallContactObject.getAsJsonPrimitive("name").getAsString())
                .zip(!overallContactObject.has("zip") ? null : overallContactObject.getAsJsonPrimitive("zip").getAsString())
                .url(!overallContactObject.has("url") ? null : overallContactObject.getAsJsonPrimitive("url").getAsString())
                .street(!overallContactObject.has("street") ? null : overallContactObject.getAsJsonPrimitive("street").getAsString())
                .type(!overallContactObject.has("type") ? null : overallContactObject.getAsJsonPrimitive("type").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchTrialsGeo createMolecularMatchTrialsGeo(@NotNull JsonObject geoObject) {
        Set<String> keysGeo = geoObject.getAsJsonObject().keySet();
        if (!EXPECTED_MOLECULARMATCH_TRAILS_GEO_ELEMENT_SIZES.contains(keysGeo.size())) {
            LOGGER.warn("Found {} in molecular match trials geo rather than the expected {}",
                    keysGeo.size(),
                    EXPECTED_MOLECULARMATCH_TRAILS_GEO_ELEMENT_SIZES);
            LOGGER.warn(keysGeo);
        }

        return ImmutableMolecularMatchTrialsGeo.builder()
                .lat(geoObject.getAsJsonPrimitive("lat").getAsString())
                .lon(geoObject.getAsJsonPrimitive("lon").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchTrialsLocation createMolecularMatchTrialsLocation(@NotNull JsonObject locationObject) {
        Set<String> keysLocation = locationObject.getAsJsonObject().keySet();
        if (!EXPECTED_MOLECULARMATCH_TRAILS_LOCATION_ELEMENT_SIZES.contains(keysLocation.size())) {
            LOGGER.warn("Found {} in molecular match trials location rather than the expected {}",
                    keysLocation.size(),
                    EXPECTED_MOLECULARMATCH_TRAILS_LOCATION_ELEMENT_SIZES);
            LOGGER.warn(keysLocation);
        }

        return ImmutableMolecularMatchTrialsLocation.builder()
                .type(locationObject.getAsJsonPrimitive("type").getAsString())
                .coordinates(jsonArrayToStringList(locationObject.getAsJsonArray("coordinates")))
                .build();
    }

    @NotNull
    private static MolecularMatchTrialsContact createMolecularMatchTrialsContact(@NotNull JsonObject contactObject) {
        Set<String> keysContact = contactObject.getAsJsonObject().keySet();
        if (!EXPECTED_MOLECULARMATCH_TRAILS_CONTACT_ELEMENT_SIZES.contains(keysContact.size())) {
            LOGGER.warn("Found {} in molecular match trials contact rather than the expected {}",
                    keysContact.size(),
                    EXPECTED_MOLECULARMATCH_TRAILS_CONTACT_ELEMENT_SIZES);
            LOGGER.warn(keysContact);
        }

        return ImmutableMolecularMatchTrialsContact.builder()
                .phone(!contactObject.has("phone") ? null : contactObject.getAsJsonPrimitive("phone").getAsString())
                .name(!contactObject.has("name") ? null : contactObject.getAsJsonPrimitive("name").getAsString())
                .email(!contactObject.has("email") ? null : contactObject.getAsJsonPrimitive("email").getAsString())
                .build();
    }

    @NotNull
    private static List<MolecularMatchTrialsTags> createMolecularMatchTrialsTags(@NotNull JsonArray arrayTags) {
        List<MolecularMatchTrialsTags> tagsList = Lists.newArrayList();
        for (JsonElement tag : arrayTags) {
            Set<String> keysTags = tag.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRAILS_TAGS_ELEMENT_SIZES.contains(keysTags.size())) {
                LOGGER.warn("Found {} in molecular match trials tags rather than the expected {}",
                        keysTags.size(),
                        EXPECTED_MOLECULARMATCH_TRAILS_TAGS_ELEMENT_SIZES);
                LOGGER.warn(keysTags);
            }

            tagsList.add(ImmutableMolecularMatchTrialsTags.builder()
                    .facet(tag.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .compositeKey(tag.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .suppress(tag.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .filterType(tag.getAsJsonObject().getAsJsonPrimitive("filterType").getAsString())
                    .term(tag.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                    .custom(tag.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .priority(tag.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .alias(!tag.getAsJsonObject().has("alias") ? null : tag.getAsJsonObject().getAsJsonPrimitive("alias").getAsString())
                    .manualSuppress(!tag.getAsJsonObject().has("manualSuppress")
                            ? null
                            : tag.getAsJsonObject().getAsJsonPrimitive("manualSuppress").getAsString())
                    .generatedBy(!tag.getAsJsonObject().has("generatedBy")
                            ? null
                            : tag.getAsJsonObject().getAsJsonPrimitive("generatedBy").getAsString())
                    .generatedByTerm(!tag.getAsJsonObject().has("generatedByTerm")
                            ? null
                            : tag.getAsJsonObject().getAsJsonPrimitive("generatedByTerm").getAsString())
                    .id(!tag.getAsJsonObject().has("id") ? null : tag.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .manualPriority(!tag.getAsJsonObject().has("manualPriority")
                            ? null
                            : tag.getAsJsonObject().getAsJsonPrimitive("manualPriority").getAsString())
                    .build());
        }

        return tagsList;
    }
}
