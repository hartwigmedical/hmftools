package com.hartwig.hmftools.iclusion.api;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.IclusionTumorLocation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTrial;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTumorLocation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class IclusionApiObjectMapper {

    private IclusionApiObjectMapper() {
    }

    @NotNull
    static List<IclusionTrial> fromApiObjects(@NotNull List<IclusionObjectStudy> studies,
            @NotNull List<IclusionObjectIndication> indications, @NotNull List<IclusionObjectGene> genes,
            @NotNull List<IclusionObjectVariant> variants) {
        List<IclusionTrial> trials = Lists.newArrayList();

        for (IclusionObjectStudy study : studies) {
            trials.add(ImmutableIclusionTrial.builder()
                    .id(study.id)
                    .title(study.title)
                    .acronym(study.acronym)
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
            @Nullable Iterable<String> indicationIds) {
        if (indicationIds == null) {
            return Lists.newArrayList();
        }

        List<IclusionTumorLocation> tumorLocations = Lists.newArrayList();
        for (String id : indicationIds) {
            IclusionObjectIndication indication = findIndicationById(indications, id);

            tumorLocations.add(ImmutableIclusionTumorLocation.builder()
                    .id(indication.id)
                    .parentId(indication.parentId)
                    .doid(indication.doid)
                    .doid2(indication.doid2)
                    .indicationName(indication.indicationName)
                    .indicationNameFull(indication.indicationNameFull)
                    .nodeIds(indication.nodeIds)
                    .build());
        }

        return tumorLocations;
    }

    @NotNull
    private static IclusionObjectIndication findIndicationById(@NotNull Iterable<IclusionObjectIndication> indications,
            @NotNull String id) {
        for (IclusionObjectIndication indication : indications) {
            if (indication.id.equals(id)) {
                return indication;
            }
        }
        throw new IllegalStateException("Could not find iClusion indication with ID '" + id + "'");
    }

    @NotNull
    private static List<IclusionMutation> buildMutations(@NotNull Iterable<IclusionObjectGene> genes,
            @NotNull Iterable<IclusionObjectVariant> variants, @NotNull Iterable<IclusionObjectMutation> mutationObjects) {
        List<IclusionMutation> mutations = Lists.newArrayList();
        for (IclusionObjectMutation mutationObject : mutationObjects) {
            IclusionObjectGene gene = findGeneById(genes, mutationObject.geneId);
            IclusionObjectVariant variant = findVariantById(variants, mutationObject.variantId);

            mutations.add(ImmutableIclusionMutation.builder().name(variant.variantName).gene(gene.geneName).build());
        }

        return mutations;
    }

    @NotNull
    private static IclusionObjectGene findGeneById(@NotNull Iterable<IclusionObjectGene> genes, @NotNull String id) {
        for (IclusionObjectGene gene : genes) {
            if (gene.id.equals(id)) {
                return gene;
            }
        }
        throw new IllegalStateException("Could not find iClusion gene with ID '" + id + "'");
    }

    @NotNull
    private static IclusionObjectVariant findVariantById(@NotNull Iterable<IclusionObjectVariant> variants, @NotNull String id) {
        for (IclusionObjectVariant variant : variants) {
            if (variant.id.equals(id)) {
                return variant;
            }
        }
        throw new IllegalStateException("Could not find iClusion variant with ID '" + id + "'");
    }
}
