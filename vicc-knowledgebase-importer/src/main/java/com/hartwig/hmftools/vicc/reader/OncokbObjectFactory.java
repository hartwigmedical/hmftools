package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.common.utils.json.JsonFunctions.nullableString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.string;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.stringList;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKb;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbClinical;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbDrugAbstract;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbGene;
import com.hartwig.hmftools.vicc.datamodel.oncokb.ImmutableOncoKbVariant;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKb;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbClinical;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbDrugAbstract;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbGene;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbVariant;

import org.jetbrains.annotations.NotNull;

final class OncokbObjectFactory {

    private OncokbObjectFactory() {
    }

    @NotNull
    static OncoKb create(@NotNull JsonObject oncoKbObject) {
        ViccDatamodelCheckerFactory.oncoKbEntryChecker().check(oncoKbObject);

        if (oncoKbObject.has("biological")) {
            return createOncoKbBiological(oncoKbObject.getAsJsonObject("biological"));
        } else if (oncoKbObject.has("clinical")) {
            return createOncoKbClinical(oncoKbObject.getAsJsonObject("clinical"));
        } else {
            throw new IllegalStateException("OncoKb object neither biological nor clinical: " + oncoKbObject);
        }
    }

    @NotNull
    private static OncoKb createOncoKbBiological(@NotNull JsonObject oncoKbBiologicalObject) {
        return ImmutableOncoKb.builder().oncoKbBiological(createBiological(oncoKbBiologicalObject)).build();
    }

    @NotNull
    private static OncoKb createOncoKbClinical(@NotNull JsonObject oncoKbClinicalObject) {
        return ImmutableOncoKb.builder().oncoKbClinical(createClinical(oncoKbClinicalObject)).build();
    }

    @NotNull
    private static OncoKbClinical createClinical(@NotNull JsonObject clinicalObject) {
        ViccDatamodelCheckerFactory.oncoKbClinicalChecker().check(clinicalObject);

        return ImmutableOncoKbClinical.builder()
                .gene(string(clinicalObject, "gene"))
                .entrezGeneId(string(clinicalObject, "Entrez Gene ID"))
                .isoform(string(clinicalObject, "Isoform"))
                .refSeq(string(clinicalObject, "RefSeq"))
                .variant(createVariant(clinicalObject.getAsJsonObject("variant")))
                .cancerType(string(clinicalObject, "cancerType"))
                .drug(string(clinicalObject, "drug"))
                .drugPmids(string(clinicalObject, "drugPmids"))
                .drugAbstracts(createDrugsAbstracts(clinicalObject.getAsJsonArray("drugAbstracts")))
                .level(string(clinicalObject, "level"))
                .levelLabel(string(clinicalObject, "level_label"))
                .build();
    }

    @NotNull
    private static List<OncoKbDrugAbstract> createDrugsAbstracts(@NotNull JsonArray drugAbstractArray) {
        List<OncoKbDrugAbstract> drugAbstractList = Lists.newArrayList();
        JsonDatamodelChecker drugAbstractChecker = ViccDatamodelCheckerFactory.oncoKbDrugsAbstractChecker();

        for (JsonElement drugAbstractElement : drugAbstractArray) {
            JsonObject drugAbstractObject = drugAbstractElement.getAsJsonObject();
            drugAbstractChecker.check(drugAbstractObject);

            drugAbstractList.add(ImmutableOncoKbDrugAbstract.builder()
                    .text(string(drugAbstractObject, "text"))
                    .link(string(drugAbstractObject, "link"))
                    .build());
        }

        return drugAbstractList;
    }

    @NotNull
    private static OncoKbBiological createBiological(@NotNull JsonObject biologicalObject) {
        ViccDatamodelCheckerFactory.oncoKbBiologicalChecker().check(biologicalObject);

        return ImmutableOncoKbBiological.builder()
                .gene(string(biologicalObject, "gene"))
                .entrezGeneId(string(biologicalObject, "Entrez Gene ID"))
                .isoform(string(biologicalObject, "Isoform"))
                .refSeq(string(biologicalObject, "RefSeq"))
                .oncokbVariant(createVariant(biologicalObject.getAsJsonObject("variant")))
                .oncogenic(string(biologicalObject, "oncogenic"))
                .mutationEffect(string(biologicalObject, "mutationEffect"))
                .mutationEffectPmids(string(biologicalObject, "mutationEffectPmids"))
                .mutationEffectAbstracts(string(biologicalObject, "mutationEffectAbstracts"))
                .build();
    }

    @NotNull
    private static OncoKbVariant createVariant(@NotNull JsonObject variantObject) {
        ViccDatamodelCheckerFactory.oncoKbVariantChecker().check(variantObject);

        return ImmutableOncoKbVariant.builder()
                .name(string(variantObject, "name"))
                .alteration(string(variantObject, "alteration"))
                .consequence(createConsequence(variantObject.getAsJsonObject("consequence")))
                .gene(createGene(variantObject.getAsJsonObject("gene")))
                .proteinStart(string(variantObject, "proteinStart"))
                .proteinEnd(string(variantObject, "proteinEnd"))
                .refResidues(nullableString(variantObject, "refResidues"))
                .variantResidues(nullableString(variantObject, "variantResidues"))
                .build();
    }

    @NotNull
    private static OncoKbConsequence createConsequence(@NotNull JsonObject consequenceObject) {
        ViccDatamodelCheckerFactory.oncoKbConsequenceChecker().check(consequenceObject);

        return ImmutableOncoKbConsequence.builder()
                .term(string(consequenceObject, "term"))
                .description(string(consequenceObject, "description"))
                .isGenerallyTruncating(string(consequenceObject, "isGenerallyTruncating"))
                .build();
    }

    @NotNull
    private static OncoKbGene createGene(@NotNull JsonObject geneObject) {
        ViccDatamodelCheckerFactory.oncoKbGeneChecker().check(geneObject);

        return ImmutableOncoKbGene.builder()
                .hugoSymbol(string(geneObject, "hugoSymbol"))
                .geneAliases(stringList(geneObject, "geneAliases"))
                .name(string(geneObject, "name"))
                .entrezGeneId(string(geneObject, "entrezGeneId"))
                .curatedIsoform(nullableString(geneObject, "curatedIsoform"))
                .curatedRefSeq(nullableString(geneObject, "curatedRefSeq"))
                .oncogene(string(geneObject, "oncogene"))
                .tsg(string(geneObject, "tsg"))
                .build();
    }
}
