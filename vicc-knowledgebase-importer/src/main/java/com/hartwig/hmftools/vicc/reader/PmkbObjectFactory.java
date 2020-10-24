package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.common.utils.json.JsonFunctions.nullableString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.string;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.pmkb.ImmutablePmkb;
import com.hartwig.hmftools.vicc.datamodel.pmkb.ImmutablePmkbGene;
import com.hartwig.hmftools.vicc.datamodel.pmkb.ImmutablePmkbTissue;
import com.hartwig.hmftools.vicc.datamodel.pmkb.ImmutablePmkbTumor;
import com.hartwig.hmftools.vicc.datamodel.pmkb.ImmutablePmkbVariant;
import com.hartwig.hmftools.vicc.datamodel.pmkb.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.pmkb.PmkbGene;
import com.hartwig.hmftools.vicc.datamodel.pmkb.PmkbTissue;
import com.hartwig.hmftools.vicc.datamodel.pmkb.PmkbTumor;
import com.hartwig.hmftools.vicc.datamodel.pmkb.PmkbVariant;

import org.jetbrains.annotations.NotNull;

final class PmkbObjectFactory {

    private PmkbObjectFactory() {
    }

    @NotNull
    static Pmkb create(@NotNull JsonObject pmkbObject) {
        ViccDatamodelCheckerFactory.pmkbEntryChecker().check(pmkbObject);

        return ImmutablePmkb.builder()
                .tumor(createTumor(pmkbObject.getAsJsonObject("tumor")))
                .tissues(createTissues(pmkbObject.getAsJsonArray("tissues")))
                .variant(createVariant(pmkbObject.getAsJsonObject("variant")))
                .build();
    }

    @NotNull
    private static PmkbTumor createTumor(@NotNull JsonObject tumorObject) {
        ViccDatamodelCheckerFactory.pmkbTumorChecker().check(tumorObject);

        return ImmutablePmkbTumor.builder()
                .name(string(tumorObject, "name"))
                .id(string(tumorObject, "id"))
                .build();
    }

    @NotNull
    private static List<PmkbTissue> createTissues(@NotNull JsonArray tissueArray) {
        List<PmkbTissue> tissueList = Lists.newArrayList();
        ViccDatamodelChecker tissueChecker = ViccDatamodelCheckerFactory.pmkbTissueChecker();

        for (JsonElement tissueElement : tissueArray) {
            JsonObject tissueObject = tissueElement.getAsJsonObject();
            tissueChecker.check(tissueObject);

            tissueList.add(ImmutablePmkbTissue.builder()
                    .name(string(tissueObject, "name"))
                    .id(string(tissueObject, "id"))
                    .build());
        }
        return tissueList;
    }

    @NotNull
    private static PmkbVariant createVariant(@NotNull JsonObject variantObject) {
        ViccDatamodelCheckerFactory.pmkbVariantChecker().check(variantObject);

        return ImmutablePmkbVariant.builder()
                .name(nullableString(variantObject, "name"))
                .coordinates(nullableString(variantObject, "coordinates"))
                .chromosome(nullableString(variantObject, "chromosome"))
                .cytoband(nullableString(variantObject, "cytoband"))
                .gene(createGene(variantObject.getAsJsonObject("gene")))
                .transcript(string(variantObject, "transcript"))
                .effect(nullableString(variantObject, "effect"))
                .codons(nullableString(variantObject, "codons"))
                .exons(nullableString(variantObject, "exons"))
                .dnaChange(nullableString(variantObject, "dna_change"))
                .aminoAcidChange(nullableString(variantObject, "amino_acid_change"))
                .germline(nullableString(variantObject, "germline"))
                .partnerGene(nullableString(variantObject, "partner_gene"))
                .cnvType(nullableString(variantObject, "cnv_type"))
                .chromosomeBasedCnv(nullableString(variantObject, "chromosome_based_cnv"))
                .variantType(nullableString(variantObject, "variant_type"))
                .cosmic(nullableString(variantObject, "cosmic"))
                .description(nullableString(variantObject, "description"))
                .descriptionType(nullableString(variantObject, "description_type"))
                .notes(nullableString(variantObject, "notes"))
                .id(nullableString(variantObject, "id"))
                .build();
    }

    @NotNull
    private static PmkbGene createGene(@NotNull JsonObject geneObject) {
        ViccDatamodelCheckerFactory.pmkbGeneChecker().check(geneObject);

        return ImmutablePmkbGene.builder()
                .name(string(geneObject, "name"))
                .createdAt(string(geneObject, "created_at"))
                .updatedAt(string(geneObject, "updated_at"))
                .activeInd(string(geneObject, "active_ind"))
                .description(nullableString(geneObject, "description"))
                .externalId(string(geneObject, "external_id"))
                .id(string(geneObject, "id"))
                .build();
    }
}
