package com.hartwig.hmftools.iclusion.api;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.IclusionTumorLocation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTrial;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTumorLocation;
import com.hartwig.hmftools.iclusion.io.IclusionTrialFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class IclusionApiObjectMapper {

    private static final Logger LOGGER = LogManager.getLogger(IclusionApiObjectMapper.class);

    private IclusionApiObjectMapper() {
    }

    @NotNull
    static List<IclusionTrial> fromApiObjects(@NotNull List<IclusionObjectStudy> studies,
            @NotNull List<IclusionObjectIndication> indications, @NotNull List<IclusionObjectGene> genes,
            @NotNull List<IclusionObjectVariant> variants) {
        List<IclusionTrial> trials = Lists.newArrayList();

        for (IclusionObjectStudy study : studies) {
            // We loose MAIN_FIELD_DELIMITER in titles for good but that is considered acceptable compromise.
            if (study.title.contains(IclusionTrialFile.MAIN_FIELD_DELIMITER)) {
                LOGGER.info("Removing trial file delimiter '{}' from study with acronym '{}'",
                        IclusionTrialFile.MAIN_FIELD_DELIMITER,
                        study.acronym);
            }

            trials.add(ImmutableIclusionTrial.builder()
                    .id(study.id)
                    .acronym(study.acronym)
                    .title(study.title.replace(IclusionTrialFile.MAIN_FIELD_DELIMITER, " "))
                    .eudra(study.eudra)
                    .nct(study.nct)
                    .ipn(study.ipn)
                    .ccmo(study.ccmo)
                    .tumorLocations(buildTumorLocations(indications, study.indicationIds))
                    .mutations(buildMutations(genes, variants, study.mutations))
                    .build());
        }

        return trials;
    }

    @NotNull
    private static List<IclusionTumorLocation> buildTumorLocations(@NotNull Iterable<IclusionObjectIndication> indications,
            @NotNull Iterable<String> indicationIds) {
        List<IclusionTumorLocation> tumorLocations = Lists.newArrayList();
        for (String id : indicationIds) {
            IclusionObjectIndication indication = findIndicationById(indications, id);
            if (indication != null) {
                List<String> doids = Lists.newArrayList();
                if (indication.doid != null) {
                    doids.add(indication.doid);
                }
                if (indication.doid2 != null) {
                    doids.add(indication.doid2);
                }

                tumorLocations.add(ImmutableIclusionTumorLocation.builder()
                        .primaryTumorLocation(indication.indicationNameFull)
                        .doids(doids)
                        .build());
            } else {
                LOGGER.warn("Could not find indication with ID '{}' in list of indications!", id);
            }
        }

        return tumorLocations;
    }

    @Nullable
    private static IclusionObjectIndication findIndicationById(@NotNull Iterable<IclusionObjectIndication> indications,
            @NotNull String id) {
        for (IclusionObjectIndication indication : indications) {
            if (indication.id.equals(id)) {
                return indication;
            }
        }
        return null;
    }

    @NotNull
    private static List<IclusionMutation> buildMutations(@NotNull Iterable<IclusionObjectGene> genes,
            @NotNull Iterable<IclusionObjectVariant> variants, @NotNull Iterable<IclusionObjectMutation> mutationObjects) {
        List<IclusionMutation> mutations = Lists.newArrayList();
        for (IclusionObjectMutation mutationObject : mutationObjects) {
            IclusionObjectGene gene = findGeneById(genes, mutationObject.geneId);
            if (gene == null) {
                LOGGER.warn("Could not find gene with ID '{}' in list of genes!", mutationObject.geneId);
            }

            IclusionObjectVariant variant = findVariantById(variants, mutationObject.variantId);
            if (variant == null) {
                LOGGER.warn("Could not find variant with ID '{}' in list of variants!", mutationObject.variantId);
            }

            if (gene != null && variant != null) {
                mutations.add(ImmutableIclusionMutation.builder().name(variant.variantName).gene(gene.geneName).build());
            }
        }

        return mutations;
    }

    @Nullable
    private static IclusionObjectGene findGeneById(@NotNull Iterable<IclusionObjectGene> genes, @NotNull String id) {
        for (IclusionObjectGene gene : genes) {
            if (gene.id.equals(id)) {
                return gene;
            }
        }
        return null;
    }

    @Nullable
    private static IclusionObjectVariant findVariantById(@NotNull Iterable<IclusionObjectVariant> variants, @NotNull String id) {
        for (IclusionObjectVariant variant : variants) {
            if (variant.id.equals(id)) {
                return variant;
            }
        }
        return null;
    }
}
