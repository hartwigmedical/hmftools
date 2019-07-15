package com.hartwig.hmftools.patientreporter;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;
import com.hartwig.hmftools.patientreporter.genepanel.ImmutableGeneModel;
import com.hartwig.hmftools.patientreporter.variants.ImmutableReportableVariant;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModelTestFactory;
import com.hartwig.hmftools.patientreporter.variants.germline.ImmutableGermlineVariant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PatientReporterTestFactory {

    private static final String ONCOGENE = "ONCO";
    private static final String TSG = "TSG";

    private PatientReporterTestFactory() {
    }

    @NotNull
    public static ImmutableGeneCopyNumber.Builder createTestCopyNumberBuilder() {
        return ImmutableGeneCopyNumber.builder()
                .start(1)
                .end(2)
                .gene(Strings.EMPTY)
                .chromosome("1")
                .chromosomeBand(Strings.EMPTY)
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(1)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .somaticRegions(1)
                .minCopyNumber(0)
                .maxCopyNumber(0)
                .transcriptID("trans")
                .transcriptVersion(0)
                .minMinorAllelePloidy(0);
    }

    @NotNull
    public static ImmutableEnrichedSomaticVariant.Builder createTestEnrichedSomaticVariantBuilder() {
        return ImmutableEnrichedSomaticVariant.builder()
                .trinucleotideContext(Strings.EMPTY)
                .highConfidenceRegion(false)
                .microhomology(Strings.EMPTY)
                .repeatSequence(Strings.EMPTY)
                .repeatCount(0)
                .clonality(Clonality.UNKNOWN)
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
                .genesEffected(0)
                .worstEffect(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.NONE)
                .worstEffectTranscript(Strings.EMPTY)
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .hotspot(Hotspot.NON_HOTSPOT)
                .recovered(false)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAllelePloidy(0)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .ploidy(0)
                .biallelic(true)
                .mappability(0);
    }

    @NotNull
    public static ImmutableGermlineVariant.Builder createTestGermlineVariantBuilder() {
        return ImmutableGermlineVariant.builder()
                .passFilter(true)
                .gene(Strings.EMPTY)
                .chromosome("1")
                .position(600)
                .hgvsCodingImpact(Strings.EMPTY)
                .hgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAllelePloidy(1D)
                .biallelic(false);
    }

    @NotNull
    public static ImmutableReportableVariant.Builder createTestReportableVariantBuilder() {
        return ImmutableReportableVariant.builder()
                .gene(Strings.EMPTY)
                .isDrupActionable(false)
                .hgvsCodingImpact(Strings.EMPTY)
                .hgvsProteinImpact(Strings.EMPTY)
                .hotspot(Hotspot.HOTSPOT)
                .clonality(Clonality.CLONAL)
                .alleleReadCount(0)
                .totalReadCount(0)
                .adjustedCopyNumber(0)
                .minorAllelePloidy(0)
                .biallelic(false)
                .adjustedVAF(0D)
                .driverCategory(DriverCategory.ONCO)
                .driverLikelihood(0D)
                .notifyClinicalGeneticist(false)
                .gDNA("1:600");
    }

    @NotNull
    public static GermlineReportingModel createTestGermlineGenesReporting() {
        Map<String, Boolean> germlineGenesReportingMap = Maps.newHashMap();
        germlineGenesReportingMap.put(ONCOGENE, true);
        germlineGenesReportingMap.put(TSG, false);
        return GermlineReportingModelTestFactory.buildFromMap(germlineGenesReportingMap);
    }

    @NotNull
    public static GermlineReportingModel createTestEmptyGermlineGenesReporting() {
        Map<String, Boolean> germlineGenesReportingMap = Maps.newHashMap();
        return GermlineReportingModelTestFactory.buildFromMap(germlineGenesReportingMap);
    }

    @NotNull
    public static GeneModel createTestGeneModel(@NotNull String oncogene, @NotNull String tsg) {
        return ImmutableGeneModel.builder()
                .addOncoDriverGenes(oncogene)
                .addTsgDriverGenes(tsg)
                .build();
    }

    @NotNull
    public static DriverCatalog createTestDriverCatalogEntry(@NotNull String gene) {
        return ImmutableDriverCatalog.builder()
                .gene(gene)
                .category(DriverCategory.ONCO)
                .driver(DriverType.MUTATION)
                .likelihoodMethod(LikelihoodMethod.NONE)
                .driverLikelihood(0D)
                .dndsLikelihood(0D)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .biallelic(false)
                .build();
    }
}
