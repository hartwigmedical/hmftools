package com.hartwig.hmftools.vicc.reader;

import java.util.List;
import java.util.Set;

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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class PmkbObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(PmkbObjectFactory.class);

    private static final List<Integer> EXPECTED_PMKB_ELEMENT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_PMKB_TUMOR_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_PMKB_TISSUE_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_PMKB_VARIANT_ELEMENT_SIZES = Lists.newArrayList(21);
    private static final List<Integer> EXPECTED_PMKB_GENE_ELEMENT_SIZES = Lists.newArrayList(7);

    private PmkbObjectFactory() {
    }

    @NotNull
    static Pmkb create(@NotNull JsonObject objectPmkb) {
        Set<String> keysPmkb = objectPmkb.keySet();
        if (!EXPECTED_PMKB_ELEMENT_SIZES.contains(keysPmkb.size())) {
            LOGGER.warn("Found {} in pmkb rather than the expected {}", keysPmkb.size(), EXPECTED_PMKB_ELEMENT_SIZES);
            LOGGER.warn(keysPmkb);
        }

        JsonObject tumor = objectPmkb.getAsJsonObject("tumor");
        JsonArray tissue = objectPmkb.getAsJsonArray("tissues");
        JsonObject variant = objectPmkb.getAsJsonObject("variant");

        return ImmutablePmkb.builder().tumor(createTumor(tumor)).tissue(createTissue(tissue)).variant(createVariantPmkb(variant)).build();
    }

    @NotNull
    private static List<PmkbTumor> createTumor(@NotNull JsonObject tumor) {
        Set<String> keysTumor = tumor.keySet();
        if (!EXPECTED_PMKB_TUMOR_ELEMENT_SIZES.contains(keysTumor.size())) {
            LOGGER.warn("Found {} in pmkb tumor rather than the expected {}", keysTumor.size(), EXPECTED_PMKB_TUMOR_ELEMENT_SIZES);
            LOGGER.warn(keysTumor);
        }

        List<PmkbTumor> listTumor = Lists.newArrayList();
        listTumor.add(ImmutablePmkbTumor.builder()
                .id(tumor.getAsJsonPrimitive("id").getAsString())
                .name(tumor.getAsJsonPrimitive("name").getAsString())
                .build());

        return listTumor;
    }

    @NotNull
    private static List<PmkbTissue> createTissue(@NotNull JsonArray tissues) {
        List<PmkbTissue> listTissue = Lists.newArrayList();
        for (JsonElement tissue : tissues) {
            Set<String> keysTissue = tissue.getAsJsonObject().keySet();
            if (!EXPECTED_PMKB_TISSUE_ELEMENT_SIZES.contains(keysTissue.size())) {
                LOGGER.warn("Found {} in pmkb tissue rather than the expected {}", keysTissue.size(), EXPECTED_PMKB_TISSUE_ELEMENT_SIZES);
                LOGGER.warn(keysTissue);
            }
            listTissue.add(ImmutablePmkbTissue.builder()
                    .id(tissue.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(tissue.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return listTissue;
    }

    @NotNull
    private static List<PmkbVariant> createVariantPmkb(@NotNull JsonObject variant) {
        Set<String> keysVariant = variant.keySet();
        if (!EXPECTED_PMKB_VARIANT_ELEMENT_SIZES.contains(keysVariant.size())) {
            LOGGER.warn("Found {} in pmkb variant rather than the expected {}", keysVariant.size(), EXPECTED_PMKB_VARIANT_ELEMENT_SIZES);
            LOGGER.warn(keysVariant);
        }

        return Lists.newArrayList(ImmutablePmkbVariant.builder()
                .aminoAcidChange(variant.get("amino_acid_change").isJsonNull() ? null : variant.get("amino_acid_change").getAsString())
                .germline(variant.get("germline").isJsonNull() ? null : variant.get("germline").getAsString())
                .partnerGene(variant.get("partner_gene").isJsonNull() ? null : variant.get("partner_gene").getAsString())
                .codons(variant.get("codons").isJsonNull() ? null : variant.get("codons").getAsString())
                .description(variant.get("description").isJsonNull() ? null : variant.get("description").getAsString())
                .exons(variant.get("exons").isJsonNull() ? null : variant.get("exons").getAsString())
                .notes(variant.get("notes").isJsonNull() ? null : variant.get("notes").getAsString())
                .cosmic(variant.get("cosmic").isJsonNull() ? null : variant.get("cosmic").getAsString())
                .effect(variant.get("effect").isJsonNull() ? null : variant.get("effect").getAsString())
                .cnvType(variant.get("cnv_type").isJsonNull() ? null : variant.get("cnv_type").getAsString())
                .id(variant.get("id").isJsonNull() ? null : variant.get("id").getAsString())
                .cytoband(variant.get("cytoband").isJsonNull() ? null : variant.get("cytoband").getAsString())
                .variantType(variant.get("variant_type").isJsonNull() ? null : variant.get("variant_type").getAsString())
                .dnaChange(variant.get("dna_change").isJsonNull() ? null : variant.get("dna_change").getAsString())
                .coordinates(variant.get("coordinates").isJsonNull() ? null : variant.get("coordinates").getAsString())
                .chromosomeBasedCnv(variant.get("chromosome_based_cnv").isJsonNull()
                        ? null
                        : variant.get("chromosome_based_cnv").getAsString())
                .gene(createGene(variant))
                .transcript(variant.getAsJsonPrimitive("transcript").getAsString())
                .descriptionType(variant.get("description_type").isJsonNull() ? null : variant.get("description_type").getAsString())
                .chromosome(variant.get("chromosome").isJsonNull() ? null : variant.get("chromosome").getAsString())
                .name(variant.get("name").isJsonNull() ? null : variant.get("name").getAsString())
                .build());
    }

    @NotNull
    private static List<PmkbGene> createGene(@NotNull JsonObject variant) {
        JsonObject gene = variant.getAsJsonObject("gene");

        Set<String> keysGene = gene.keySet();
        if (!EXPECTED_PMKB_GENE_ELEMENT_SIZES.contains(keysGene.size())) {
            LOGGER.warn("Found {} in pmkb gene rather than the expected {}", keysGene.size(), EXPECTED_PMKB_GENE_ELEMENT_SIZES);
            LOGGER.warn(keysGene);
        }

        List<PmkbGene> listGene = Lists.newArrayList();
        listGene.add(ImmutablePmkbGene.builder()
                .description(gene.get("description").isJsonNull() ? null : gene.getAsJsonPrimitive("description").getAsString())
                .createdAt(gene.getAsJsonPrimitive("created_at").getAsString())
                .updatedAt(gene.getAsJsonPrimitive("updated_at").getAsString())
                .activeInd(gene.getAsJsonPrimitive("active_ind").getAsString())
                .externalId(gene.getAsJsonPrimitive("external_id").getAsString())
                .id(gene.getAsJsonPrimitive("id").getAsString())
                .name(gene.getAsJsonPrimitive("name").getAsString())
                .build());
        return listGene;
    }
}
