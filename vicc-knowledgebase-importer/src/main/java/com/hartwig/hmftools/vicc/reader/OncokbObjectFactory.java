package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.nullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.string;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.stringList;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbClinical;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbDrugAbstracts;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncokb;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncokbGene;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncokbVariant;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbClinical;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbDrugAbstracts;
import com.hartwig.hmftools.vicc.datamodel.oncokb.Oncokb;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncokbGene;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncokbVariant;

import org.jetbrains.annotations.NotNull;

final class OncokbObjectFactory {

    private OncokbObjectFactory() {
    }

    @NotNull
    static Oncokb create(@NotNull JsonObject oncoKbObject) {
        ViccDatamodelCheckerFactory.oncoKbEntryChecker().check(oncoKbObject);

        if (oncoKbObject.has("biological")) {
            return createOncoKbBiological(oncoKbObject);
        } else if (oncoKbObject.has("clinical")) {
            return createOncoKbClinical(oncoKbObject);
        } else {
            throw new IllegalStateException("OncoKb object neither biological nor clinical: " + oncoKbObject);
        }
    }

    @NotNull
    private static Oncokb createOncoKbBiological(@NotNull JsonObject oncoKbObject) {
        return ImmutableOncokb.builder().oncoKbBiological(createBiologicalOncoKb(oncoKbObject.getAsJsonObject("biological"))).build();
    }

    @NotNull
    private static Oncokb createOncoKbClinical(@NotNull JsonObject oncoKbObject) {
        return ImmutableOncokb.builder().oncoKbClinical(createClinicalOncoKb(oncoKbObject.getAsJsonObject("clinical"))).build();
    }

    @NotNull
    private static OncoKbClinical createClinicalOncoKb(@NotNull JsonObject clinicalObject) {
        ViccDatamodelCheckerFactory.oncoKbClinicalChecker().check(clinicalObject);

        return ImmutableOncoKbClinical.builder()
                .refSeq(string(clinicalObject, "RefSeq"))
                .level(string(clinicalObject, "level"))
                .isoform(string(clinicalObject, "Isoform"))
                .oncokbVariant(createVariantOncoKb(clinicalObject.getAsJsonObject("variant")))
                .entrezGeneID(string(clinicalObject, "Entrez Gene ID"))
                .drugPmids(string(clinicalObject, "drugPmids"))
                .cancerType(string(clinicalObject, "cancerType"))
                .drug(string(clinicalObject, "drug"))
                .gene(string(clinicalObject, "gene"))
                .levelLabel(string(clinicalObject, "level_label"))
                .oncoKbDrugAbstracts(createDrugsAbstracts(clinicalObject.getAsJsonArray("drugAbstracts")))
                .build();
    }

    @NotNull
    private static List<OncoKbDrugAbstracts> createDrugsAbstracts(@NotNull JsonArray drugsAbstractArray) {
        List<OncoKbDrugAbstracts> drugsAbstractList = Lists.newArrayList();
        for (JsonElement drugAbstracts : drugsAbstractArray) {
            JsonObject drugAbstractsObject = drugAbstracts.getAsJsonObject();

            ViccDatamodelCheckerFactory.oncoKbDrugsAbstractChecker().check(drugAbstractsObject);
            drugsAbstractList.add(ImmutableOncoKbDrugAbstracts.builder()
                    .text(string(drugAbstractsObject, "text"))
                    .link(string(drugAbstractsObject, "link"))
                    .build());
        }

        return drugsAbstractList;
    }

    @NotNull
    private static OncoKbBiological createBiologicalOncoKb(@NotNull JsonObject biologicalObject) {
        ViccDatamodelCheckerFactory.oncoKbBiologicalChecker().check(biologicalObject);

        return ImmutableOncoKbBiological.builder()
                .mutationEffectPmids(string(biologicalObject, "mutationEffectPmids"))
                .isoform(string(biologicalObject, "Isoform"))
                .oncokbVariant(createVariantOncoKb(biologicalObject.getAsJsonObject("variant")))
                .entrezGeneID(string(biologicalObject, "Entrez Gene ID"))
                .oncogenic(string(biologicalObject, "oncogenic"))
                .mutationEffect(string(biologicalObject, "mutationEffect"))
                .refSeq(string(biologicalObject, "RefSeq"))
                .gene(string(biologicalObject, "gene"))
                .mutationEffectAbstracts(string(biologicalObject, "mutationEffectAbstracts"))
                .build();
    }

    @NotNull
    private static OncokbVariant createVariantOncoKb(@NotNull JsonObject variantObject) {
        ViccDatamodelCheckerFactory.oncoKbVariantChecker().check(variantObject);

        return ImmutableOncokbVariant.builder()
                .variantResidues(nullableString(variantObject, "variantResidues"))
                .proteinStart(string(variantObject, "proteinStart"))
                .name(string(variantObject, "name"))
                .proteinEnd(string(variantObject, "proteinEnd"))
                .refResidues(nullableString(variantObject, "refResidues"))
                .alteration(string(variantObject, "alteration"))
                .oncoKbConsequence(createConsequenceOncokb(variantObject.getAsJsonObject("consequence")))
                .oncokbGene(createGeneOncoKb(variantObject.getAsJsonObject("gene")))
                .build();
    }

    @NotNull
    private static OncoKbConsequence createConsequenceOncokb(@NotNull JsonObject consequenceObject) {
        ViccDatamodelCheckerFactory.oncoKbConsequenceChecker().check(consequenceObject);

        return ImmutableOncoKbConsequence.builder()
                .term(string(consequenceObject, "term"))
                .description(string(consequenceObject, "description"))
                .isGenerallyTruncating(string(consequenceObject, "isGenerallyTruncating"))
                .build();
    }

    @NotNull
    private static OncokbGene createGeneOncoKb(@NotNull JsonObject geneObject) {
        ViccDatamodelCheckerFactory.oncoKbGeneChecker().check(geneObject);

        return ImmutableOncokbGene.builder()
                .oncogene(string(geneObject, "oncogene"))
                .name(string(geneObject, "name"))
                .hugoSymbol(string(geneObject, "hugoSymbol"))
                .curatedRefSeq(nullableString(geneObject, "curatedRefSeq"))
                .entrezGeneId(string(geneObject, "entrezGeneId"))
                .geneAliases(stringList(geneObject, "geneAliases"))
                .tsg(string(geneObject, "tsg"))
                .curatedIsoform(nullableString(geneObject, "curatedIsoform"))
                .build();
    }
}
