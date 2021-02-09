package com.hartwig.hmftools.ckb.reader.indication;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.datamodel.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.datamodel.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableClinicalTrialInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableEvidenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.common.IndicationInfo;
import com.hartwig.hmftools.ckb.datamodel.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.indication.ImmutableIndication;
import com.hartwig.hmftools.ckb.datamodel.indication.Indication;
import com.hartwig.hmftools.ckb.reader.CkbJsonDirectoryReader;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class IndicationReader extends CkbJsonDirectoryReader<Indication> {

    public IndicationReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected Indication read(@NotNull final JsonObject object) {
        JsonDatamodelChecker indicationChecker = IndicationDataModelChecker.indicationObjectChecker();
        indicationChecker.check(object);

        return ImmutableIndication.builder()
                .id(JsonFunctions.integer(object, "id"))
                .name(JsonFunctions.string(object, "name"))
                .source(JsonFunctions.string(object, "source"))
                .definition(JsonFunctions.nullableString(object, "definition"))
                .currentPreferredTerm(JsonFunctions.nullableString(object, "currentPreferredTerm"))
                .lastUpdateDateFromDO(JsonFunctions.nullableString(object, "lastUpdateDateFromDO"))
                .altId(JsonFunctions.stringList(object, "altIds"))
                .termId(JsonFunctions.string(object, "termId"))
                .evidence(extractEvidence(object.getAsJsonArray("evidence")))
                .clinicalTrial(extractClinicalTrials(object.getAsJsonArray("clinicalTrials")))
                .build();
    }

    @NotNull
    private static List<EvidenceInfo> extractEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> evidences = Lists.newArrayList();
        JsonDatamodelChecker indicationEvidenceChecker = IndicationDataModelChecker.indicationEvidenceObjectChecker();

        for (JsonElement evidence : jsonArray) {
            JsonObject evidenceJsonObject = evidence.getAsJsonObject();
            indicationEvidenceChecker.check(evidenceJsonObject);

            evidences.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(evidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(evidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(evidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(evidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(evidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(evidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(evidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(evidenceJsonObject, "responseType"))
                    .reference(extractReference(evidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidenceJsonObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationEvidenceMolecularProfileChecker = IndicationDataModelChecker.indicationEvidenceMolecularProfileObjectChecker();
        indicationEvidenceMolecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationEvidenceTherapyChecker = IndicationDataModelChecker.indicationEvidenceTherapyObjectChecker();
        indicationEvidenceTherapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationEvidenceIndicationChecker = IndicationDataModelChecker.indicationEvidenceIndicationObjectChecker();
        indicationEvidenceIndicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static List<ReferenceInfo> extractReference(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker indicationEvidenceReferenceChecker = IndicationDataModelChecker.indicationEvidenceReferenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();
            indicationEvidenceReferenceChecker.check(referenceJsonObject);

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
        JsonDatamodelChecker indicationClinicalTrialChecker = IndicationDataModelChecker.indicationClinicaltrialObjectChecker();

        for (JsonElement clinicalTrial : jsonArray) {
            JsonObject clinicalTrialJsonObject = clinicalTrial.getAsJsonObject();
            indicationClinicalTrialChecker.check(clinicalTrialJsonObject);

            clinicalTrials.add(ImmutableClinicalTrialInfo.builder()
                    .nctId(JsonFunctions.string(clinicalTrialJsonObject, "nctId"))
                    .title(JsonFunctions.string(clinicalTrialJsonObject, "title"))
                    .phase(JsonFunctions.string(clinicalTrialJsonObject, "phase"))
                    .recruitment(JsonFunctions.string(clinicalTrialJsonObject, "recruitment"))
                    .therapy(extractTherapyList(clinicalTrialJsonObject.getAsJsonArray("therapies")))
                    .build());
        }
        return clinicalTrials;
    }

    @NotNull
    private static List<TherapyInfo> extractTherapyList(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker indicationClinicalTrialTherapiesChecker = IndicationDataModelChecker.indicationEvidenceTherapyObjectChecker();

        for (JsonElement therapy : jsonArray) {
            JsonObject therapyJsonObject = therapy.getAsJsonObject();
            indicationClinicalTrialTherapiesChecker.check(therapyJsonObject);

            therapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(therapyJsonObject, "id"))
                    .therapyName(JsonFunctions.string(therapyJsonObject, "therapyName"))
                    .synonyms(JsonFunctions.nullableString(therapyJsonObject, "synonyms"))
                    .build());
        }
        return therapies;
    }
}
