package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.jsonArrayToStringList;

import java.util.List;
import java.util.Set;

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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class OncokbObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(OncokbObjectFactory.class);

    private static final List<Integer> EXPECTED_ONCOKB_ELEMENT_SIZES = Lists.newArrayList(1);
    private static final List<Integer> EXPECTED_ONCOKB_CLINICAL_ELEMENT_SIZES = Lists.newArrayList(11);
    private static final List<Integer> EXPECTED_ONCOKB_DRUGS_ABSTRACT_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_ONCOKB_BIOLOGICAL_ELEMENT_SIZES = Lists.newArrayList(9);
    private static final List<Integer> EXPECTED_ONCOKB_VARIANT_ELEMENT_SIZES = Lists.newArrayList(8);
    private static final List<Integer> EXPECTED_ONCOKB_CONSEQUENCE_ELEMENT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_ONCOKB_GENE_ELEMENT_SIZES = Lists.newArrayList(8);

    private OncokbObjectFactory() {
    }

    @NotNull
    static Oncokb create(@NotNull JsonObject objectOncoKb) {
        Set<String> keysOncokb = objectOncoKb.keySet();
        if (!EXPECTED_ONCOKB_ELEMENT_SIZES.contains(keysOncokb.size())) {
            LOGGER.warn("Found {} in oncokb rather than the expected {}", keysOncokb.size(), EXPECTED_ONCOKB_ELEMENT_SIZES);
            LOGGER.warn(keysOncokb);
        }

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
        Set<String> keysClinical = objectClinical.keySet();
        if (!EXPECTED_ONCOKB_CLINICAL_ELEMENT_SIZES.contains(keysClinical.size())) {
            LOGGER.warn("Found {} in oncokb clinical rather than the expected {}",
                    keysClinical.size(),
                    EXPECTED_ONCOKB_CLINICAL_ELEMENT_SIZES);
            LOGGER.warn(keysClinical);
        }

        return ImmutableOncoKbClinical.builder()
                .RefSeq(objectClinical.getAsJsonPrimitive("RefSeq").getAsString())
                .level(objectClinical.getAsJsonPrimitive("level").getAsString())
                .Isoform(objectClinical.getAsJsonPrimitive("Isoform").getAsString())
                .oncokbVariant(createVariantOncoKb(objectClinical.getAsJsonObject("variant")))
                .entrezGeneID(objectClinical.getAsJsonPrimitive("Entrez Gene ID").getAsString())
                .drugPmids(objectClinical.getAsJsonPrimitive("drugPmids").getAsString())
                .cancerType(objectClinical.getAsJsonPrimitive("cancerType").getAsString())
                .drug(objectClinical.getAsJsonPrimitive("drug").getAsString())
                .gene(objectClinical.getAsJsonPrimitive("gene").getAsString())
                .levelLabel(objectClinical.getAsJsonPrimitive("level_label").getAsString())
                .oncoKbDrugAbstracts(createDrugsAbstracts(objectClinical.getAsJsonArray("drugAbstracts")))
                .build();
    }

    @NotNull
    private static List<OncoKbDrugAbstracts> createDrugsAbstracts(@NotNull JsonArray arrayDrugsAbstracts) {
        List<OncoKbDrugAbstracts> listDrugsAbstracts = Lists.newArrayList();
        for (JsonElement drugAbstracts : arrayDrugsAbstracts) {
            Set<String> keysBiological = drugAbstracts.getAsJsonObject().keySet();

            if (!EXPECTED_ONCOKB_DRUGS_ABSTRACT_ELEMENT_SIZES.contains(keysBiological.size())) {
                LOGGER.warn("Found {} in oncokb drugs abstracts rather than the expected {}",
                        keysBiological.size(),
                        EXPECTED_ONCOKB_DRUGS_ABSTRACT_ELEMENT_SIZES);
                LOGGER.warn(keysBiological);
            }
            listDrugsAbstracts.add(ImmutableOncoKbDrugAbstracts.builder()
                    .text(drugAbstracts.getAsJsonObject().getAsJsonPrimitive("text").getAsString())
                    .link(drugAbstracts.getAsJsonObject().getAsJsonPrimitive("link").getAsString())
                    .build());
        }
        return listDrugsAbstracts;
    }

    @NotNull
    private static OncoKbBiological createBiologicalOncoKb(@NotNull JsonObject objectBiological) {
        Set<String> keysBiological = objectBiological.keySet();
        if (!EXPECTED_ONCOKB_BIOLOGICAL_ELEMENT_SIZES.contains(keysBiological.size())) {
            LOGGER.warn("Found {} in oncokb biological rather than the expected {}",
                    keysBiological.size(),
                    EXPECTED_ONCOKB_BIOLOGICAL_ELEMENT_SIZES);
            LOGGER.warn(keysBiological);
        }

        return ImmutableOncoKbBiological.builder()
                .mutationEffectPmids(objectBiological.getAsJsonPrimitive("mutationEffectPmids").getAsString())
                .Isoform(objectBiological.getAsJsonPrimitive("Isoform").getAsString())
                .oncokbVariant(createVariantOncoKb(objectBiological.getAsJsonObject("variant")))
                .entrezGeneID(objectBiological.getAsJsonPrimitive("Entrez Gene ID").getAsString())
                .oncogenic(objectBiological.getAsJsonPrimitive("oncogenic").getAsString())
                .mutationEffect(objectBiological.getAsJsonPrimitive("mutationEffect").getAsString())
                .RefSeq(objectBiological.getAsJsonPrimitive("RefSeq").getAsString())
                .gene(objectBiological.getAsJsonPrimitive("gene").getAsString())
                .mutationEffectAbstracts(objectBiological.getAsJsonPrimitive("mutationEffectAbstracts").getAsString())
                .build();
    }

    @NotNull
    private static OncokbVariant createVariantOncoKb(@NotNull JsonObject objectVariant) {
        Set<String> keysVariant = objectVariant.keySet();

        if (!EXPECTED_ONCOKB_VARIANT_ELEMENT_SIZES.contains(keysVariant.size())) {
            LOGGER.warn("Found {} in oncokb variant rather than the expected {}",
                    keysVariant.size(),
                    EXPECTED_ONCOKB_VARIANT_ELEMENT_SIZES);
            LOGGER.warn(keysVariant);
        }

        return ImmutableOncokbVariant.builder()
                .variantResidues(objectVariant.get("variantResidues").isJsonNull()
                        ? null
                        : objectVariant.getAsJsonPrimitive("variantResidues").getAsString())
                .proteinStart(objectVariant.getAsJsonPrimitive("proteinStart").getAsString())
                .name(objectVariant.getAsJsonPrimitive("name").getAsString())
                .proteinEnd(objectVariant.getAsJsonPrimitive("proteinEnd").getAsString())
                .refResidues(objectVariant.get("refResidues").isJsonNull()
                        ? null
                        : objectVariant.getAsJsonPrimitive("refResidues").getAsString())
                .alteration(objectVariant.getAsJsonPrimitive("alteration").getAsString())
                .oncoKbConsequence(createConsequenceOncokb(objectVariant.getAsJsonObject("consequence")))
                .oncokbGene(createGeneOncoKb(objectVariant.getAsJsonObject("gene")))
                .build();
    }

    @NotNull
    private static OncoKbConsequence createConsequenceOncokb(@NotNull JsonObject objectConsequence) {
        Set<String> keysConsequence = objectConsequence.keySet();

        if (!EXPECTED_ONCOKB_CONSEQUENCE_ELEMENT_SIZES.contains(keysConsequence.size())) {
            LOGGER.warn("Found {} in oncokb consequence rather than the expected {}",
                    keysConsequence.size(),
                    EXPECTED_ONCOKB_CONSEQUENCE_ELEMENT_SIZES);
            LOGGER.warn(keysConsequence);
        }

        return ImmutableOncoKbConsequence.builder()
                .term(objectConsequence.getAsJsonPrimitive("term").getAsString())
                .description(objectConsequence.getAsJsonPrimitive("description").getAsString())
                .isGenerallyTruncating(objectConsequence.getAsJsonPrimitive("isGenerallyTruncating").getAsString())
                .build();
    }

    @NotNull
    private static OncokbGene createGeneOncoKb(@NotNull JsonObject objectGene) {
        Set<String> keysGene = objectGene.keySet();

        if (!EXPECTED_ONCOKB_GENE_ELEMENT_SIZES.contains(keysGene.size())) {
            LOGGER.warn("Found {} in oncokb gene rather than the expected {}", keysGene.size(), EXPECTED_ONCOKB_GENE_ELEMENT_SIZES);
            LOGGER.warn(keysGene);
        }

        return ImmutableOncokbGene.builder()
                .oncogene(objectGene.getAsJsonPrimitive("oncogene").getAsString())
                .name(objectGene.getAsJsonPrimitive("name").getAsString())
                .hugoSymbol(objectGene.getAsJsonPrimitive("hugoSymbol").getAsString())
                .curatedRefSeq(objectGene.get("curatedRefSeq").isJsonNull()
                        ? null
                        : objectGene.getAsJsonPrimitive("curatedRefSeq").getAsString())
                .entrezGeneId(objectGene.getAsJsonPrimitive("entrezGeneId").getAsString())
                .geneAliases(Lists.newArrayList(jsonArrayToStringList(objectGene.getAsJsonArray("geneAliases"))))
                .tsg(objectGene.getAsJsonPrimitive("tsg").getAsString())
                .curatedIsoform(objectGene.get("curatedIsoform").isJsonNull()
                        ? null
                        : objectGene.getAsJsonPrimitive("curatedIsoform").getAsString())
                .build();
    }
}
