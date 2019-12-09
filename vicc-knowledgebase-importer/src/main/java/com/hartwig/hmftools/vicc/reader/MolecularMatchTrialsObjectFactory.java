package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.nullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalJsonObject;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalNullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalStringList;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.string;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.stringList;

import java.util.List;

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

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class MolecularMatchTrialsObjectFactory {

    private MolecularMatchTrialsObjectFactory() {
    }

    @NotNull
    static MolecularMatchTrials create(@NotNull JsonObject molecularMatchTrialsObject) {
        ViccDatamodelCheckerFactory.molecularMatchTrialsEntryChecker().check(molecularMatchTrialsObject);

        return ImmutableMolecularMatchTrials.builder()
                .status(string(molecularMatchTrialsObject, "status"))
                .startDate(optionalNullableString(molecularMatchTrialsObject, "startDate"))
                .title(string(molecularMatchTrialsObject, "title"))
                .molecularAlterations(stringList(molecularMatchTrialsObject, "molecularAlterations"))
                .score(string(molecularMatchTrialsObject, "_score"))
                .intervation(createInterventions(molecularMatchTrialsObject.getAsJsonArray("interventions")))
                .locations(createLocations(molecularMatchTrialsObject.getAsJsonArray("locations")))
                .briefTitle(nullableString(molecularMatchTrialsObject, "briefTitle"))
                .overallContact(createOverallContact(optionalJsonObject(molecularMatchTrialsObject, "overallContact")))
                .link(string(molecularMatchTrialsObject, "link"))
                .phase(string(molecularMatchTrialsObject, "phase"))
                .tags(createTags(molecularMatchTrialsObject.getAsJsonArray("tags")))
                .id(string(molecularMatchTrialsObject, "id"))
                .studyType(string(molecularMatchTrialsObject, "studyType"))
                .build();
    }

    @NotNull
    private static List<MolecularMatchTrialsIntervation> createInterventions(@NotNull JsonArray interventionArray) {
        List<MolecularMatchTrialsIntervation> molecularMatchTrialsInterventionList = Lists.newArrayList();
        ViccDatamodelChecker interventionChecker = ViccDatamodelCheckerFactory.molecularMatchTrialsInterventionChecker();

        for (JsonElement interventionElement : interventionArray) {
            JsonObject interventionObject = interventionElement.getAsJsonObject();
            interventionChecker.check(interventionObject);

            molecularMatchTrialsInterventionList.add(ImmutableMolecularMatchTrialsIntervation.builder()
                    .intervention_name(optionalString(interventionObject, "intervention_name"))
                    .other_name(optionalStringList(interventionObject, "other_name"))
                    .description(optionalString(interventionObject, "description"))
                    .arm_group_label(optionalStringList(interventionObject, "arm_group_label"))
                    .intervention_type(optionalString(interventionObject, "intervention_type"))
                    .build());
        }

        return molecularMatchTrialsInterventionList;
    }

    @NotNull
    private static List<MolecularMatchTrialsLocations> createLocations(@NotNull JsonArray locationArray) {
        List<MolecularMatchTrialsLocations> locationList = Lists.newArrayList();
        ViccDatamodelChecker locationChecker = ViccDatamodelCheckerFactory.molecularMatchTrialsLocationChecker();

        for (JsonElement locationElement : locationArray) {
            JsonObject locationObject = locationElement.getAsJsonObject();
            locationChecker.check(locationObject);

            locationList.add(ImmutableMolecularMatchTrialsLocations.builder()
                    .status(string(locationObject, "status"))
                    .last_name(optionalString(locationObject, "last_name"))
                    .email(optionalString(locationObject, "email"))
                    .phone(optionalString(locationObject, "phone"))
                    .phone_backup(optionalString(locationObject, "phone_backup"))
                    .email_backup(optionalString(locationObject, "email_backup"))
                    .last_name_backup(optionalString(locationObject, "last_name_backup"))
                    .phone_ext_backup(optionalString(locationObject, "phone_ext_backup"))
                    .phone_ext(optionalString(locationObject, "phone_ext"))
                    .city(optionalString(locationObject, "city"))
                    .valid(optionalString(locationObject, "_valid"))
                    .zip(optionalString(locationObject, "zip"))
                    .created(optionalString(locationObject, "created"))
                    .country(optionalString(locationObject, "country"))
                    .number(optionalString(locationObject, "number"))
                    .id(optionalString(locationObject, "id"))
                    .lastUpdated(optionalString(locationObject, "lastUpdated"))
                    .contact(createContact(optionalJsonObject(locationObject, "contact")))
                    .state(optionalString(locationObject, "state"))
                    .street(optionalString(locationObject, "street"))
                    .location(createSubLocation(optionalJsonObject(locationObject, "location")))
                    .po_box(optionalString(locationObject, "po_box"))
                    .failedGeocode(optionalString(locationObject, "failedGeocode"))
                    .geo(createGeo(optionalJsonObject(locationObject, "geo")))
                    .validMessage(optionalString(locationObject, "_validMessage"))
                    .name(optionalString(locationObject, "name"))
                    .build());
        }

        return locationList;
    }

    @Nullable
    private static MolecularMatchTrialsOverallContact createOverallContact(@Nullable JsonObject overallContactObject) {
        if (overallContactObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchTrialsOverallContactChecker().check(overallContactObject);

        return ImmutableMolecularMatchTrialsOverallContact.builder()
                .phone(optionalNullableString(overallContactObject, "phone"))
                .last_name(optionalString(overallContactObject, "last_name"))
                .email(optionalString(overallContactObject, "email"))
                .affiliation(optionalNullableString(overallContactObject, "affiliation"))
                .phone_ext(optionalString(overallContactObject, "phone_ext"))
                .country(optionalString(overallContactObject, "country"))
                .city(optionalString(overallContactObject, "city"))
                .name(optionalString(overallContactObject, "name"))
                .zip(optionalString(overallContactObject, "zip"))
                .url(optionalString(overallContactObject, "url"))
                .street(optionalString(overallContactObject, "street"))
                .type(optionalString(overallContactObject, "type"))
                .build();
    }

    @Nullable
    private static MolecularMatchTrialsGeo createGeo(@Nullable JsonObject geoObject) {
        if (geoObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchTrialsGeoChecker().check(geoObject);

        return ImmutableMolecularMatchTrialsGeo.builder().lat(string(geoObject, "lat")).lon(string(geoObject, "lon")).build();
    }

    @Nullable
    private static MolecularMatchTrialsLocation createSubLocation(@Nullable JsonObject subLocationObject) {
        if (subLocationObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchTrialsSubLocationChecker().check(subLocationObject);

        return ImmutableMolecularMatchTrialsLocation.builder()
                .type(string(subLocationObject, "type"))
                .coordinates(stringList(subLocationObject, "coordinates"))
                .build();
    }

    @Nullable
    private static MolecularMatchTrialsContact createContact(@Nullable JsonObject contactObject) {
        if (contactObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchTrialsContactChecker().check(contactObject);

        return ImmutableMolecularMatchTrialsContact.builder()
                .phone(optionalString(contactObject, "phone"))
                .name(optionalString(contactObject, "name"))
                .email(optionalString(contactObject, "email"))
                .build();
    }

    @NotNull
    private static List<MolecularMatchTrialsTags> createTags(@NotNull JsonArray tagArray) {
        List<MolecularMatchTrialsTags> tagList = Lists.newArrayList();
        ViccDatamodelChecker tagChecker = ViccDatamodelCheckerFactory.molecularMatchTrialsTagChecker();

        for (JsonElement tagElement : tagArray) {
            JsonObject tagObject = tagElement.getAsJsonObject();
            tagChecker.check(tagObject);

            tagList.add(ImmutableMolecularMatchTrialsTags.builder()
                    .facet(string(tagObject, "facet"))
                    .compositeKey(string(tagObject, "compositeKey"))
                    .suppress(string(tagObject, "suppress"))
                    .filterType(string(tagObject, "filterType"))
                    .term(string(tagObject, "term"))
                    .custom(string(tagObject, "custom"))
                    .priority(string(tagObject, "priority"))
                    .alias(optionalString(tagObject, "alias"))
                    .manualSuppress(optionalString(tagObject, "manualSuppress"))
                    .generatedBy(optionalString(tagObject, "generatedBy"))
                    .generatedByTerm(optionalString(tagObject, "generatedByTerm"))
                    .id(optionalString(tagObject, "id"))
                    .manualPriority(optionalString(tagObject, "manualPriority"))
                    .build());
        }

        return tagList;
    }
}
