package com.hartwig.hmftools.common.actionability;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.PrimaryTumorToDOIDMapping;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberEvidenceAnalyzerFactory;
import com.hartwig.hmftools.common.actionability.fusion.FusionEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.fusion.FusionEvidenceAnalyzerFactory;
import com.hartwig.hmftools.common.actionability.somaticvariant.SomaticVariantEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.somaticvariant.SomaticVariantEvidenceAnalyzerFactory;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ActionabilityAnalyzer {

    private static final String ACTIONABLE_VARIANT_FILE = "actionableVariants.tsv";
    private static final String ACTIONABLE_RANGES_FILE = "actionableRanges.tsv";
    private static final String ACTIONABLE_CNV_FILE = "actionableCNVs.tsv";
    private static final String ACTIONABLE_FUSION_PAIR_FILE = "actionableFusionPairs.tsv";
    private static final String ACTIONABLE_PROMISCUOUS_FIVE_FILE = "actionablePromiscuousFive.tsv";
    private static final String ACTIONABLE_PROMISCUOUS_THREE_FILE = "actionablePromiscuousThree.tsv";

    private static final String CANCER_TYPE_DOID_MAPPING_FILE = "knowledgebaseCancerTypes.tsv";

    @NotNull
    private final SomaticVariantEvidenceAnalyzer variantAnalyzer;
    @NotNull
    private final CopyNumberEvidenceAnalyzer cnvAnalyzer;
    @NotNull
    private final FusionEvidenceAnalyzer fusionAnalyzer;
    @NotNull
    private final CancerTypeAnalyzer cancerTypeAnalyzer;
    @NotNull
    private final PrimaryTumorToDOIDMapping primaryTumorToDOIDMapping;

    @NotNull
    public static ActionabilityAnalyzer fromKnowledgebase(@NotNull String knowledgebasePath) throws IOException {
        String basePath = knowledgebasePath + File.separator;
        SomaticVariantEvidenceAnalyzer variantAnalyzer =
                SomaticVariantEvidenceAnalyzerFactory.loadFromFileVariantsAndFileRanges(basePath + ACTIONABLE_VARIANT_FILE,
                        basePath + ACTIONABLE_RANGES_FILE);

        CopyNumberEvidenceAnalyzer cnvAnalyzer = CopyNumberEvidenceAnalyzerFactory.loadFromFileCNVs(basePath + ACTIONABLE_CNV_FILE);

        FusionEvidenceAnalyzer fusionAnalyzer = FusionEvidenceAnalyzerFactory.loadFromFileFusions(basePath + ACTIONABLE_FUSION_PAIR_FILE,
                basePath + ACTIONABLE_PROMISCUOUS_FIVE_FILE,
                basePath + ACTIONABLE_PROMISCUOUS_THREE_FILE);

        CancerTypeAnalyzer cancerTypeAnalyzer =
                CancerTypeAnalyzer.loadFromFile(knowledgebasePath + File.separator + CANCER_TYPE_DOID_MAPPING_FILE);

        return new ActionabilityAnalyzer(variantAnalyzer,
                cnvAnalyzer,
                fusionAnalyzer,
                cancerTypeAnalyzer,
                PrimaryTumorToDOIDMapping.createFromResource());
    }

    private ActionabilityAnalyzer(@NotNull final SomaticVariantEvidenceAnalyzer variantAnalyzer,
            @NotNull final CopyNumberEvidenceAnalyzer cnvAnalyzer, @NotNull final FusionEvidenceAnalyzer fusionAnalyzer,
            @NotNull final CancerTypeAnalyzer cancerTypeAnalyzer, @NotNull final PrimaryTumorToDOIDMapping primaryTumorToDOIDMapping) {
        this.variantAnalyzer = variantAnalyzer;
        this.cnvAnalyzer = cnvAnalyzer;
        this.fusionAnalyzer = fusionAnalyzer;
        this.cancerTypeAnalyzer = cancerTypeAnalyzer;
        this.primaryTumorToDOIDMapping = primaryTumorToDOIDMapping;
    }

    @VisibleForTesting
    @NotNull
    SomaticVariantEvidenceAnalyzer variantAnalyzer() {
        return variantAnalyzer;
    }

    @VisibleForTesting
    @NotNull
    CopyNumberEvidenceAnalyzer cnvAnalyzer() {
        return cnvAnalyzer;
    }

    @VisibleForTesting
    @NotNull
    FusionEvidenceAnalyzer fusionAnalyzer() {
        return fusionAnalyzer;
    }

    @NotNull
    public Map<EnrichedSomaticVariant, List<EvidenceItem>> evidenceForSomaticVariants(@NotNull List<EnrichedSomaticVariant> variants,
            @Nullable String primaryTumorLocation) {
        Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant = Maps.newHashMap();

        List<EnrichedSomaticVariant> variantsOnActionableGenes = variants.stream()
                .filter(variant -> variantAnalyzer.actionableGenes().contains(variant.gene()))
                .collect(Collectors.toList());

        String doidsPrimaryTumorLocation = doidForPrimaryTumorLocation(primaryTumorLocation);

        for (EnrichedSomaticVariant actionableVariant : variantsOnActionableGenes) {
            evidencePerVariant.put(actionableVariant,
                    variantAnalyzer.evidenceForSomaticVariant(actionableVariant, doidsPrimaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerVariant;
    }

    @NotNull
    public Map<GeneCopyNumber, List<EvidenceItem>> evidenceForCopyNumbers(@NotNull List<GeneCopyNumber> geneCopyNumbers,
            @Nullable String primaryTumorLocation) {
        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerCopyNumber = Maps.newHashMap();

        // TODO (KODU): Should also filter on significant event rather than assume caller has filtered already.
        List<GeneCopyNumber> geneCopyNumbersOnActionableGenes = geneCopyNumbers.stream()
                .filter(geneCopyNumber -> cnvAnalyzer.actionableGenes().contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        String doidsPrimaryTumorLocation = doidForPrimaryTumorLocation(primaryTumorLocation);

        for (GeneCopyNumber geneCopyNumber : geneCopyNumbersOnActionableGenes) {
            evidencePerCopyNumber.put(geneCopyNumber,
                    cnvAnalyzer.evidenceForCopyNumber(geneCopyNumber, doidsPrimaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerCopyNumber;
    }

    @NotNull
    public Map<GeneFusion, List<EvidenceItem>> evidenceForFusions(@NotNull List<GeneFusion> fusions,
            @Nullable String primaryTumorLocation) {
        Map<GeneFusion, List<EvidenceItem>> evidencePerFusion = Maps.newHashMap();

        List<GeneFusion> fusionsOnActionableGenes = fusions.stream()
                .filter(fusion -> fusionAnalyzer.actionableGenes().contains(fusion.upstreamLinkedAnnotation().geneName())
                        || fusionAnalyzer.actionableGenes().contains(fusion.downstreamLinkedAnnotation().geneName()))
                .collect(Collectors.toList());

        String doidsPrimaryTumorLocation = doidForPrimaryTumorLocation(primaryTumorLocation);

        // TODO (KODU): Should reuse "favor canonical" rules from SV analyser here but have re-implemented for now.
        for (GeneFusion actionableFusion : uniqueGeneFusions(fusionsOnActionableGenes)) {
            evidencePerFusion.put(actionableFusion,
                    fusionAnalyzer.evidenceForFusion(actionableFusion, doidsPrimaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerFusion;
    }

    @NotNull
    private static List<GeneFusion> uniqueGeneFusions(@NotNull List<GeneFusion> fusions) {
        List<GeneFusion> allUniqueGeneFusions = Lists.newArrayList();
        Map<FiveThreePair, List<GeneFusion>> fusionsPerFiveThreePair = Maps.newHashMap();
        for (GeneFusion fusion : fusions) {
            FiveThreePair key =
                    new FiveThreePair(fusion.upstreamLinkedAnnotation().geneName(), fusion.downstreamLinkedAnnotation().geneName());
            List<GeneFusion> fusionsForKey = fusionsPerFiveThreePair.get(key);
            if (fusionsForKey == null) {
                fusionsForKey = Lists.newArrayList();
            }
            fusionsForKey.add(fusion);
            fusionsPerFiveThreePair.put(key, fusionsForKey);
        }

        for (Map.Entry<FiveThreePair, List<GeneFusion>> fusionsPerKey : fusionsPerFiveThreePair.entrySet()) {
            allUniqueGeneFusions.add(favorCanonical(fusionsPerKey.getValue()));
        }

        return allUniqueGeneFusions;
    }

    @NotNull
    private static GeneFusion favorCanonical(@NotNull List<GeneFusion> fusions) {
        for (GeneFusion fusion : fusions) {
            if (fusion.upstreamLinkedAnnotation().isCanonical() && fusion.downstreamLinkedAnnotation().isCanonical()) {
                return fusion;
            }
        }

        // KODU: If there is no canonical-canonical fusion, return the first one arbitrarily.
        assert !fusions.isEmpty();
        return fusions.get(0);
    }

    @Nullable
    private String doidForPrimaryTumorLocation(@Nullable String primaryTumorLocation) {
        if (primaryTumorLocation == null) {
            return null;
        }

        Set<String> doidsPrimaryTumorLocation = primaryTumorToDOIDMapping.doidsForPrimaryTumorLocation(primaryTumorLocation);

        // TODO (KODU): Deal with sets of DOIDs. Also every non-null tumor location should return a DOID.
        return doidsPrimaryTumorLocation != null && !doidsPrimaryTumorLocation.isEmpty()
                ? doidsPrimaryTumorLocation.iterator().next()
                : null;
    }

    private static class FiveThreePair {

        @NotNull
        private final String fiveGene;
        @NotNull
        private final String threeGene;

        private FiveThreePair(@NotNull final String fiveGene, @NotNull final String threeGene) {
            this.fiveGene = fiveGene;
            this.threeGene = threeGene;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final FiveThreePair that = (FiveThreePair) o;
            return Objects.equals(fiveGene, that.fiveGene) && Objects.equals(threeGene, that.threeGene);
        }

        @Override
        public int hashCode() {
            return Objects.hash(fiveGene, threeGene);
        }
    }
}
