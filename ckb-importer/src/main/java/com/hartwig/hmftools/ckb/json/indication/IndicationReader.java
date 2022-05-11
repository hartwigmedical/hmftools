package com.hartwig.hmftools.ckb.json.indication;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.json.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableEvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableReferenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.util.DateConverter;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class IndicationReader extends CkbJsonDirectoryReader<JsonIndication> {
    private static final Logger LOGGER = LogManager.getLogger(IndicationReader.class);

    public IndicationReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected JsonIndication read(@NotNull final JsonObject object) {
        JsonDatamodelChecker indicationChecker = IndicationDataModelChecker.indicationObjectChecker();
        indicationChecker.check(object);

        return ImmutableJsonIndication.builder()
                .id(JsonFunctions.integer(object, "id"))
                .name(JsonFunctions.string(object, "name"))
                .source(JsonFunctions.string(object, "source"))
                .definition(JsonFunctions.nullableString(object, "definition"))
                .currentPreferredTerm(JsonFunctions.nullableString(object, "currentPreferredTerm"))
                .lastUpdateDateFromDO(DateConverter.toDate(JsonFunctions.nullableString(object, "lastUpdateDateFromDO")))
                .altIds(JsonFunctions.stringList(object, "altIds"))
                .termId(JsonFunctions.string(object, "termId"))
                .evidence(extractEvidence(object.getAsJsonArray("evidence")))
                .clinicalTrials(extractClinicalTrials(object.getAsJsonArray("clinicalTrials")))
                .build();
    }

    @NotNull
    private static List<EvidenceInfo> extractEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> evidences = Lists.newArrayList();
        JsonDatamodelChecker evidenceChecker = IndicationDataModelChecker.evidenceObjectChecker();

        for (JsonElement evidence : jsonArray) {
            JsonObject evidenceJsonObject = evidence.getAsJsonObject();
            evidenceChecker.check(evidenceJsonObject);

            evidences.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(evidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(evidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(evidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(evidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(evidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(evidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(evidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(evidenceJsonObject, "responseType"))
                    .references(extractReferences(evidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidenceJsonObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker evidenceMolecularProfileChecker = IndicationDataModelChecker.evidenceMolecularProfileObjectChecker();
        evidenceMolecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker evidenceTherapyChecker = IndicationDataModelChecker.evidenceTherapyObjectChecker();
        evidenceTherapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.optionalStringList(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker evidenceIndicationChecker = IndicationDataModelChecker.evidenceIndicationObjectChecker();
        evidenceIndicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static List<ReferenceInfo> extractReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker evidenceReferenceChecker = IndicationDataModelChecker.evidenceReferenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();
            evidenceReferenceChecker.check(referenceJsonObject);

            references.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.integer(referenceJsonObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceJsonObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(referenceJsonObject, "title"))
                    .url(JsonFunctions.nullableString(referenceJsonObject, "url"))
                    .build());
        }
        return references;
    }

    @NotNull
    private static List<ClinicalTrialInfo> extractClinicalTrials(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialInfo> clinicalTrials = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrialChecker = IndicationDataModelChecker.clinicalTrialObjectChecker();

        for (JsonElement clinicalTrial : jsonArray) {
            JsonObject clinicalTrialJsonObject = clinicalTrial.getAsJsonObject();
            clinicalTrialChecker.check(clinicalTrialJsonObject);

            String nctId = JsonFunctions.string(clinicalTrialJsonObject, "nctId");
            String phase = JsonFunctions.nullableString(clinicalTrialJsonObject, "phase");

            if (phase == null) {
                LOGGER.warn("phase of study '{}' is null in IndicationReader", nctId);
            }

            clinicalTrials.add(ImmutableClinicalTrialInfo.builder()
                    .nctId(nctId)
                    .title(JsonFunctions.string(clinicalTrialJsonObject, "title"))
                    .phase(phase)
                    .recruitment(JsonFunctions.string(clinicalTrialJsonObject, "recruitment"))
                    .therapies(extractTherapies(clinicalTrialJsonObject.getAsJsonArray("therapies")))
                    .build());
        }
        return clinicalTrials;
    }

    @NotNull
    private static List<TherapyInfo> extractTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrialTherapyChecker = IndicationDataModelChecker.evidenceTherapyObjectChecker();

        for (JsonElement therapy : jsonArray) {
            JsonObject therapyJsonObject = therapy.getAsJsonObject();
            clinicalTrialTherapyChecker.check(therapyJsonObject);

            therapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(therapyJsonObject, "id"))
                    .therapyName(JsonFunctions.string(therapyJsonObject, "therapyName"))
                    .synonyms(JsonFunctions.optionalStringList(therapyJsonObject, "synonyms"))
                    .build());
        }
        return therapies;
    }
}