package com.hartwig.hmftools.protect;

import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.germline.ImmutableReportableGermlineVariant;
import com.hartwig.hmftools.protect.homozygousdisruption.ImmutableReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.structural.ImmutableReportableGeneDisruption;
import com.hartwig.hmftools.protect.variants.germline.ConditionReportingVariant;
import com.hartwig.hmftools.protect.variants.germline.GermlineReporting;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.variants.germline.ImmutableGermlineReporting;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ProtectTestFactory {

    private static final String KNOWLEDGEBASE_DIRECTORY = Resources.getResource("actionability").getPath();

    public static final String ONCOGENE = "KIT";
    public static final String TSG = "PTEN";

    private ProtectTestFactory() {
    }

    @NotNull
    public static ActionabilityAnalyzer loadTestActionabilityAnalyzer() {
        try {
            return ActionabilityAnalyzer.fromKnowledgebase(KNOWLEDGEBASE_DIRECTORY);
        } catch (IOException e) {
            throw new IllegalStateException(e);
        }
    }

    @NotNull
    public static ImmutableSomaticVariantImpl.Builder createTestSomaticVariantBuilder() {
        return ImmutableSomaticVariantImpl.builder()
                .qual(100)
                .trinucleotideContext(Strings.EMPTY)
                .microhomology(Strings.EMPTY)
                .repeatSequence(Strings.EMPTY)
                .repeatCount(0)
                .kataegis(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .type(VariantType.UNDEFINED)
                .filter("PASS")
                .totalReadCount(0)
                .alleleReadCount(0)
                .gene(Strings.EMPTY)
                .genesAffected(0)
                .worstEffect(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.NONE)
                .worstEffectTranscript(Strings.EMPTY)
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .hotspot(Hotspot.NON_HOTSPOT)
                .recovered(false)
                .reported(false)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAlleleCopyNumber(0)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .variantCopyNumber(0)
                .biallelic(true)
                .subclonalLikelihood(0)
                .tier(VariantTier.UNKNOWN)
                .mappability(0);
    }

    @NotNull
    public static ImmutableReportableGermlineVariant.Builder createTestGermlineVariantBuilder() {
        return ImmutableReportableGermlineVariant.builder()
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .codingEffect(CodingEffect.UNDEFINED)
                .hgvsCoding(Strings.EMPTY)
                .alleleReadCount(0)
                .totalReadCount(0)
                .biallelic(false);
    }

    @NotNull
    public static ImmutableReportableGainLoss.Builder createTestReportableGainLossBuilder() {
        return ImmutableReportableGainLoss.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .copies(0)
                .interpretation(CopyNumberInterpretation.FULL_LOSS);
    }

    @NotNull
    public static ImmutableReportableHomozygousDisruption.Builder createTestReportableHomozygousDisruptionBuilder() {
        return ImmutableReportableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(Strings.EMPTY);
    }

    @NotNull
    public static ImmutableReportableGeneDisruption.Builder createTestReportableGeneDisruptionBuilder() {
        return ImmutableReportableGeneDisruption.builder()
                .gene(Strings.EMPTY)
                .location(Strings.EMPTY)
                .range(Strings.EMPTY)
                .type(Strings.EMPTY)
                .firstAffectedExon(0)
                .junctionCopyNumber(0D)
                .undisruptedCopyNumber(0D);
    }

    @NotNull
    public static ImmutableEvidenceItem.Builder createTestEvidenceBuilder() {
        return ImmutableEvidenceItem.builder()
                .event(Strings.EMPTY)
                .source(ActionabilitySource.CIVIC)
                .reference(Strings.EMPTY)
                .drug(Strings.EMPTY)
                .drugsType(Strings.EMPTY)
                .level(EvidenceLevel.LEVEL_A)
                .response(Strings.EMPTY)
                .isOnLabel(false)
                .cancerType(Strings.EMPTY)
                .scope(EvidenceScope.SPECIFIC);
    }

    @NotNull
    public static GermlineReportingModel createTestGermlineGenesReportingBialleic() {
        Map<String, GermlineReporting> germlineGenesReportingMap = Maps.newHashMap();
        GermlineReporting germlineReportingTrue = ImmutableGermlineReporting.builder()
                .notifyClinicalGeneticus(true)
                .condition(ConditionReportingVariant.BIALLELIC_ONLY)
                .variant("")
                .build();
        GermlineReporting germlineReportingFalse = ImmutableGermlineReporting.builder()
                .notifyClinicalGeneticus(false)
                .condition(ConditionReportingVariant.BIALLELIC_ONLY)
                .variant("")
                .build();
        germlineGenesReportingMap.put(ONCOGENE, germlineReportingTrue);
        germlineGenesReportingMap.put(TSG, germlineReportingFalse);
        return new GermlineReportingModel(germlineGenesReportingMap);
    }

    @NotNull
    public static GermlineReportingModel createTestGermlineGenesReportingMonoalleic() {
        Map<String, GermlineReporting> germlineGenesReportingMap = Maps.newHashMap();
        GermlineReporting germlineReportingTrue = ImmutableGermlineReporting.builder()
                .notifyClinicalGeneticus(true)
                .condition(ConditionReportingVariant.ALL)
                .variant("")
                .build();
        GermlineReporting germlineReportingFalse = ImmutableGermlineReporting.builder()
                .notifyClinicalGeneticus(false)
                .condition(ConditionReportingVariant.ALL)
                .variant("")
                .build();
        germlineGenesReportingMap.put(ONCOGENE, germlineReportingTrue);
        germlineGenesReportingMap.put(TSG, germlineReportingFalse);
        return new GermlineReportingModel(germlineGenesReportingMap);
    }

    @NotNull
    public static GermlineReportingModel createTestEmptyGermlineGenesReporting() {
        Map<String, GermlineReporting> germlineGenesReportingMap = Maps.newHashMap();
        return new GermlineReportingModel(germlineGenesReportingMap);
    }

    @NotNull
    public static ImmutableDriverCatalog.Builder createTestDriverCatalogBuilder() {
        return ImmutableDriverCatalog.builder()
                .biallelic(false)
                .category(DriverCategory.ONCO)
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .driver(DriverType.MUTATION)
                .driverLikelihood(1.0)
                .dndsLikelihood(1.0)
                .likelihoodMethod(LikelihoodMethod.DNDS)
                .minCopyNumber(0D)
                .maxCopyNumber(0D)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0);
    }
}
