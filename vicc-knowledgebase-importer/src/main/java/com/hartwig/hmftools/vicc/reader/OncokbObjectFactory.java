package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalNullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.string;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.toStringList;

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
    static Oncokb create(@NotNull JsonObject objectOncoKb) {

        ViccDatamodelCheckerFactory.oncoKbEntryChecker().check(objectOncoKb);

        if (objectOncoKb.has("biological")) {
            return createOncoKbBiological(objectOncoKb);
        } else if (objectOncoKb.has("clinical")) {
            return createOncoKbClinical(objectOncoKb);
        } else {
            throw new IllegalStateException("OncoKb object neither biological nor clinical: " + objectOncoKb);
        }
    }

    @NotNull
    private static Oncokb createOncoKbBiological(@NotNull JsonObject objectOncoKb) {
        return ImmutableOncokb.builder().oncoKbBiological(createBiologicalOncoKb(objectOncoKb.getAsJsonObject("biological"))).build();
    }

    @NotNull
    private static Oncokb createOncoKbClinical(@NotNull JsonObject objectOncoKb) {
        return ImmutableOncokb.builder().oncoKbClinical(createClinicalOncoKb(objectOncoKb.getAsJsonObject("clinical"))).build();
    }

    @NotNull
    private static OncoKbClinical createClinicalOncoKb(@NotNull JsonObject objectClinical) {
        ViccDatamodelCheckerFactory.oncoKbClinicalChecker().check(objectClinical);

        return ImmutableOncoKbClinical.builder()
                .RefSeq(string(objectClinical, "RefSeq"))
                .level(string(objectClinical,"level"))
                .Isoform(string(objectClinical,"Isoform"))
                .oncokbVariant(createVariantOncoKb(objectClinical.getAsJsonObject("variant")))
                .entrezGeneID(string(objectClinical,"Entrez Gene ID"))
                .drugPmids(string(objectClinical,"drugPmids"))
                .cancerType(string(objectClinical,"cancerType"))
                .drug(string(objectClinical,"drug"))
                .gene(string(objectClinical,"gene"))
                .levelLabel(string(objectClinical,"level_label"))
                .oncoKbDrugAbstracts(createDrugsAbstracts(objectClinical.getAsJsonArray("drugAbstracts")))
                .build();
    }

    @NotNull
    private static List<OncoKbDrugAbstracts> createDrugsAbstracts(@NotNull JsonArray arrayDrugsAbstracts) {
        List<OncoKbDrugAbstracts> listDrugsAbstracts = Lists.newArrayList();
        for (JsonElement drugAbstracts : arrayDrugsAbstracts) {

            JsonObject drugAbstractsObject = drugAbstracts.getAsJsonObject();

            ViccDatamodelCheckerFactory.oncoKbDrugsAbstractChecker().check(drugAbstractsObject);
            listDrugsAbstracts.add(ImmutableOncoKbDrugAbstracts.builder()
                    .text(string(drugAbstractsObject, "text"))
                    .link(string(drugAbstractsObject, "link"))
                    .build());
        }
        return listDrugsAbstracts;
    }

    @NotNull
    private static OncoKbBiological createBiologicalOncoKb(@NotNull JsonObject objectBiological) {
        ViccDatamodelCheckerFactory.oncoKbBiologicalChecker().check(objectBiological);

        return ImmutableOncoKbBiological.builder()
                .mutationEffectPmids(string(objectBiological, "mutationEffectPmids"))
                .Isoform(string(objectBiological, "Isoform"))
                .oncokbVariant(createVariantOncoKb(objectBiological.getAsJsonObject("variant")))
                .entrezGeneID(string(objectBiological,"Entrez Gene ID"))
                .oncogenic(string(objectBiological,"oncogenic"))
                .mutationEffect(string(objectBiological,"mutationEffect"))
                .RefSeq(string(objectBiological,"RefSeq"))
                .gene(string(objectBiological,"gene"))
                .mutationEffectAbstracts(string(objectBiological,"mutationEffectAbstracts"))
                .build();
    }

    @NotNull
    private static OncokbVariant createVariantOncoKb(@NotNull JsonObject objectVariant) {
        ViccDatamodelCheckerFactory.oncoKbVariantChecker().check(objectVariant);

        return ImmutableOncokbVariant.builder()
                .variantResidues(optionalNullableString(objectVariant, "variantResidues"))
                .proteinStart(string(objectVariant, "proteinStart"))
                .name(string(objectVariant, "name"))
                .proteinEnd(string(objectVariant, "proteinEnd"))
                .refResidues(optionalNullableString(objectVariant, "refResidues"))
                .alteration(string(objectVariant, "alteration"))
                .oncoKbConsequence(createConsequenceOncokb(objectVariant.getAsJsonObject("consequence")))
                .oncokbGene(createGeneOncoKb(objectVariant.getAsJsonObject("gene")))
                .build();
    }

    @NotNull
    private static OncoKbConsequence createConsequenceOncokb(@NotNull JsonObject objectConsequence) {
        ViccDatamodelCheckerFactory.oncoKbConsequenceChecker().check(objectConsequence);


        return ImmutableOncoKbConsequence.builder()
                .term(string(objectConsequence, "term"))
                .description(string(objectConsequence, "description"))
                .isGenerallyTruncating(string(objectConsequence, "isGenerallyTruncating"))
                .build();
    }

    @NotNull
    private static OncokbGene createGeneOncoKb(@NotNull JsonObject objectGene) {
        ViccDatamodelCheckerFactory.oncoKbGeneChecker().check(objectGene);


        return ImmutableOncokbGene.builder()
                .oncogene(string(objectGene, "oncogene"))
                .name(string(objectGene, "name"))
                .hugoSymbol(string(objectGene, "hugoSymbol"))
                .curatedRefSeq(optionalNullableString(objectGene,"curatedRefSeq"))
                .entrezGeneId(string(objectGene, "entrezGeneId"))
                .geneAliases(toStringList(objectGene.getAsJsonArray("geneAliases")))
                .tsg(string(objectGene, "tsg"))
                .curatedIsoform(optionalNullableString(objectGene,"curatedIsoform"))
                .build();
    }
}
