package com.hartwig.hmftools.ckb.json.reference;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.json.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.json.common.DrugInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.GeneInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableDrugInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableEvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableGeneInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableReferenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableTreatmentApproachInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableVariantInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.common.TreatmentApproachInfo;
import com.hartwig.hmftools.ckb.json.common.VariantInfo;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ReferenceReader extends CkbJsonDirectoryReader<JsonReference> {

    public ReferenceReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected JsonReference read(@NotNull final JsonObject object) {
        JsonDatamodelChecker referenceChecker = ReferenceDataModelChecker.referenceObjectChecker();
        referenceChecker.check(object);

        return ImmutableJsonReference.builder()
                .id(JsonFunctions.integer(object, "id"))
                .pubMedId(JsonFunctions.nullableString(object, "pubMedId"))
                .title(JsonFunctions.nullableString(object, "title"))
                .url(JsonFunctions.nullableString(object, "url"))
                .authors(JsonFunctions.nullableString(object, "authors"))
                .journal(JsonFunctions.nullableString(object, "journal"))
                .volume(JsonFunctions.nullableString(object, "volume"))
                .issue(JsonFunctions.nullableString(object, "issue"))
                .date(JsonFunctions.nullableString(object, "date"))
                .abstractText(JsonFunctions.nullableString(object, "abstractText"))
                .year(JsonFunctions.nullableString(object, "year"))
                .drugs(extractDrugs(object.getAsJsonArray("drugs")))
                .genes(extractGenes(object.getAsJsonArray("genes")))
                .evidence(extractEvidence(object.getAsJsonArray("evidence")))
                .therapies(extractTherapies(object.getAsJsonArray("therapies")))
                .treatmentApproaches(extractTreatmentApproaches(object.getAsJsonArray("treatmentApproaches")))
                .variants(extractVariants(object.getAsJsonArray("variants")))
                .build();
    }

    @NotNull
    private static List<DrugInfo> extractDrugs(@NotNull JsonArray jsonArray) {
        List<DrugInfo> referenceDrugs = Lists.newArrayList();
        JsonDatamodelChecker drugChecker = ReferenceDataModelChecker.drugObjectChecker();

        for (JsonElement drug : jsonArray) {
            JsonObject drugJsonObject = drug.getAsJsonObject();
            drugChecker.check(drugJsonObject);

            referenceDrugs.add(ImmutableDrugInfo.builder()
                    .id(JsonFunctions.integer(drugJsonObject, "id"))
                    .drugName(JsonFunctions.string(drugJsonObject, "drugName"))
                    .terms(JsonFunctions.stringList(drugJsonObject, "terms"))
                    .build());
        }
        return referenceDrugs;
    }

    @NotNull
    private static List<GeneInfo> extractGenes(@NotNull JsonArray jsonArray) {
        List<GeneInfo> genes = Lists.newArrayList();
        JsonDatamodelChecker geneChecker = ReferenceDataModelChecker.geneObjectChecker();

        for (JsonElement gene : jsonArray) {
            JsonObject geneJsonObject = gene.getAsJsonObject();
            geneChecker.check(geneJsonObject);

            genes.add(ImmutableGeneInfo.builder()
                    .id(JsonFunctions.integer(geneJsonObject, "id"))
                    .geneSymbol(JsonFunctions.string(geneJsonObject, "geneSymbol"))
                    .terms(JsonFunctions.stringList(geneJsonObject, "terms"))
                    .build());
        }
        return genes;
    }

    @NotNull
    private static List<EvidenceInfo> extractEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> evidences = Lists.newArrayList();
        JsonDatamodelChecker evidenceChecker = ReferenceDataModelChecker.evidenceObjectChecker();

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
        JsonDatamodelChecker molecularProfileChecker = ReferenceDataModelChecker.molecularProfileObjectChecker();
        molecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyChecker = ReferenceDataModelChecker.therapyChecker();
        therapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.optionalStringList(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationChecker = ReferenceDataModelChecker.indicationChecker();
        indicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static List<ReferenceInfo> extractReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker referenceChecker = ReferenceDataModelChecker.referenceChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();
            referenceChecker.check(referenceJsonObject);

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
    private static List<TherapyInfo> extractTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker therapyChecker = ReferenceDataModelChecker.therapyChecker();

        for (JsonElement therapy : jsonArray) {
            JsonObject therapyJsonObject = therapy.getAsJsonObject();
            therapyChecker.check(therapyJsonObject);

            therapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(therapyJsonObject, "id"))
                    .therapyName(JsonFunctions.string(therapyJsonObject, "therapyName"))
                    .synonyms(JsonFunctions.optionalStringList(therapyJsonObject, "synonyms"))
                    .build());
        }
        return therapies;
    }

    @NotNull
    private static List<TreatmentApproachInfo> extractTreatmentApproaches(@NotNull JsonArray jsonArray) {
        List<TreatmentApproachInfo> treatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker treatmentApproachChecker = ReferenceDataModelChecker.treatmentApproachObjectChecker();

        for (JsonElement treatmentApproach : jsonArray) {
            JsonObject treatmentApproachJsonObject = treatmentApproach.getAsJsonObject();
            treatmentApproachChecker.check(treatmentApproachJsonObject);

            treatmentApproaches.add(ImmutableTreatmentApproachInfo.builder()
                    .id(JsonFunctions.integer(treatmentApproachJsonObject, "id"))
                    .name(JsonFunctions.string(treatmentApproachJsonObject, "name"))
                    .profileName(JsonFunctions.string(treatmentApproachJsonObject, "profileName"))
                    .build());
        }
        return treatmentApproaches;
    }

    @NotNull
    private static List<VariantInfo> extractVariants(@NotNull JsonArray jsonArray) {
        List<VariantInfo> variants = Lists.newArrayList();
        JsonDatamodelChecker variantChecker = ReferenceDataModelChecker.variantObjectChecker();

        for (JsonElement variant : jsonArray) {
            JsonObject variantJsonObject = variant.getAsJsonObject();
            variantChecker.check(variantJsonObject);

            variants.add(ImmutableVariantInfo.builder()
                    .id(JsonFunctions.integer(variantJsonObject, "id"))
                    .fullName(JsonFunctions.string(variantJsonObject, "fullName"))
                    .impact(JsonFunctions.nullableString(variantJsonObject, "impact"))
                    .proteinEffect(JsonFunctions.nullableString(variantJsonObject, "proteinEffect"))
                    .build());
        }
        return variants;
    }
}
