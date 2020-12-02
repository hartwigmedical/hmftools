package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.common.utils.json.JsonFunctions.nullableString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalJsonObject;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalNullableString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalStringList;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.string;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.stringList;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsContact;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsGeo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsIntervention;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsOverallContact;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsSubLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.ImmutableMolecularMatchTrialsTag;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsContact;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsGeo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsIntervention;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsOverallContact;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsSubLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsTag;

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
                .briefTitle(nullableString(molecularMatchTrialsObject, "briefTitle"))
                .studyType(string(molecularMatchTrialsObject, "studyType"))
                .molecularAlterations(stringList(molecularMatchTrialsObject, "molecularAlterations"))
                .score(string(molecularMatchTrialsObject, "_score"))
                .interventions(createInterventions(molecularMatchTrialsObject.getAsJsonArray("interventions")))
                .locations(createLocations(molecularMatchTrialsObject.getAsJsonArray("locations")))
                .overallContact(createOverallContact(optionalJsonObject(molecularMatchTrialsObject, "overallContact")))
                .link(string(molecularMatchTrialsObject, "link"))
                .phase(string(molecularMatchTrialsObject, "phase"))
                .tags(createTags(molecularMatchTrialsObject.getAsJsonArray("tags")))
                .id(string(molecularMatchTrialsObject, "id"))
                .build();
    }

    @NotNull
    private static List<MolecularMatchTrialsIntervention> createInterventions(@NotNull JsonArray interventionArray) {
        List<MolecularMatchTrialsIntervention> molecularMatchTrialsInterventionList = Lists.newArrayList();
        JsonDatamodelChecker interventionChecker = ViccDatamodelCheckerFactory.molecularMatchTrialsInterventionChecker();

        for (JsonElement interventionElement : interventionArray) {
            JsonObject interventionObject = interventionElement.getAsJsonObject();
            interventionChecker.check(interventionObject);

            molecularMatchTrialsInterventionList.add(ImmutableMolecularMatchTrialsIntervention.builder()
                    .interventionName(optionalString(interventionObject, "intervention_name"))
                    .otherNames(optionalStringList(interventionObject, "other_name"))
                    .interventionType(optionalString(interventionObject, "intervention_type"))
                    .armGroupLabels(optionalStringList(interventionObject, "arm_group_label"))
                    .description(optionalString(interventionObject, "description"))
                    .build());
        }

        return molecularMatchTrialsInterventionList;
    }

    @NotNull
    private static List<MolecularMatchTrialsLocation> createLocations(@NotNull JsonArray locationArray) {
        List<MolecularMatchTrialsLocation> locationList = Lists.newArrayList();
        JsonDatamodelChecker locationChecker = ViccDatamodelCheckerFactory.molecularMatchTrialsLocationChecker();

        for (JsonElement locationElement : locationArray) {
            JsonObject locationObject = locationElement.getAsJsonObject();
            locationChecker.check(locationObject);

            locationList.add(ImmutableMolecularMatchTrialsLocation.builder()
                    .status(string(locationObject, "status"))
                    .name(optionalString(locationObject, "name"))
                    .contact(createContact(optionalJsonObject(locationObject, "contact")))
                    .lastName(optionalString(locationObject, "last_name"))
                    .email(optionalString(locationObject, "email"))
                    .phone(optionalString(locationObject, "phone"))
                    .phoneExt(optionalString(locationObject, "phone_ext"))
                    .lastNameBackup(optionalString(locationObject, "last_name_backup"))
                    .emailBackup(optionalString(locationObject, "email_backup"))
                    .phoneBackup(optionalString(locationObject, "phone_backup"))
                    .phoneExtBackup(optionalString(locationObject, "phone_ext_backup"))
                    .subLocation(createSubLocation(optionalJsonObject(locationObject, "location")))
                    .street(optionalString(locationObject, "street"))
                    .city(optionalString(locationObject, "city"))
                    .zip(optionalString(locationObject, "zip"))
                    .state(optionalString(locationObject, "state"))
                    .country(optionalString(locationObject, "country"))
                    .number(optionalString(locationObject, "number"))
                    .poBox(optionalString(locationObject, "po_box"))
                    .id(optionalString(locationObject, "id"))
                    .valid(optionalString(locationObject, "_valid"))
                    .validMessage(optionalString(locationObject, "_validMessage"))
                    .created(optionalString(locationObject, "created"))
                    .lastUpdated(optionalString(locationObject, "lastUpdated"))
                    .failedGeocode(optionalString(locationObject, "failedGeocode"))
                    .geo(createGeo(optionalJsonObject(locationObject, "geo")))
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
                .name(optionalString(overallContactObject, "name"))
                .type(optionalString(overallContactObject, "type"))
                .affiliation(optionalNullableString(overallContactObject, "affiliation"))
                .lastName(optionalString(overallContactObject, "last_name"))
                .email(optionalString(overallContactObject, "email"))
                .phone(optionalNullableString(overallContactObject, "phone"))
                .phoneExt(optionalString(overallContactObject, "phone_ext"))
                .street(optionalString(overallContactObject, "street"))
                .city(optionalString(overallContactObject, "city"))
                .zip(optionalString(overallContactObject, "zip"))
                .country(optionalString(overallContactObject, "country"))
                .url(optionalString(overallContactObject, "url"))
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
    private static MolecularMatchTrialsSubLocation createSubLocation(@Nullable JsonObject subLocationObject) {
        if (subLocationObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchTrialsSubLocationChecker().check(subLocationObject);

        return ImmutableMolecularMatchTrialsSubLocation.builder()
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
                .name(optionalString(contactObject, "name"))
                .email(optionalString(contactObject, "email"))
                .phone(optionalString(contactObject, "phone"))
                .build();
    }

    @NotNull
    private static List<MolecularMatchTrialsTag> createTags(@NotNull JsonArray tagArray) {
        List<MolecularMatchTrialsTag> tagList = Lists.newArrayList();
        JsonDatamodelChecker tagChecker = ViccDatamodelCheckerFactory.molecularMatchTrialsTagChecker();

        for (JsonElement tagElement : tagArray) {
            JsonObject tagObject = tagElement.getAsJsonObject();
            tagChecker.check(tagObject);

            tagList.add(ImmutableMolecularMatchTrialsTag.builder()
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
