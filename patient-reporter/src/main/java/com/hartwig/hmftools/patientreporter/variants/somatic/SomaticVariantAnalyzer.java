package com.hartwig.hmftools.patientreporter.variants.somatic;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.OncoDrivers;
import com.hartwig.hmftools.common.drivercatalog.TsgDrivers;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SomaticVariantAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);

    private static final List<CodingEffect> TSG_CODING_EFFECTS_TO_REPORT =
            Lists.newArrayList(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE, CodingEffect.SPLICE);

    private static final Set<String> ONCO_GENES_WITH_SPLICE_EFFECTS = Sets.newHashSet("MET", "JAK2");

    private static final List<CodingEffect> ONCO_CODING_EFFECTS_TO_REPORT = Lists.newArrayList(CodingEffect.MISSENSE);

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static SomaticVariantAnalysis run(@NotNull List<EnrichedSomaticVariant> variants, @NotNull GeneModel panelGeneModel,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation) {
        double microsatelliteIndelsPerMb = MicrosatelliteAnalyzer.determineMicrosatelliteIndelsPerMb(variants);
        int tumorMutationalLoad = MutationalLoadAnalyzer.determineTumorMutationalLoad(variants);
        double tumorMutationalBurden = MutationalBurdenAnalyzer.determineTumorMutationalBurden(variants);

        List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(OncoDrivers.drivers(variants));
        driverCatalog.addAll(TsgDrivers.drivers(variants));

        List<EnrichedSomaticVariant> variantsToReport =
                variants.stream().filter(includeFilter(panelGeneModel)).collect(Collectors.toList());

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant =
                actionabilityAnalyzer.evidenceForSomaticVariants(variants, primaryTumorLocation);

        // Extract somatic evidence for high drivers variants
        Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariantHighDriver = Maps.newHashMap();
        for (Map.Entry<EnrichedSomaticVariant, List<EvidenceItem>> entry : evidencePerVariant.entrySet()) {
            String gene = entry.getKey().gene();
            for (DriverCatalog catalog : driverCatalog) {
                if (catalog.gene().equals(gene)) {
                    if (catalog.driverLikelihood() > 0.8) {
                        evidencePerVariantHighDriver.put(entry.getKey(), entry.getValue());
                    }
                }
            }
        }

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerVariantHighDriver);

        // Check that all variants with high level evidence are reported (since they are in the driver catalog).
        for (Map.Entry<EnrichedSomaticVariant, List<EvidenceItem>> entry : evidencePerVariantHighDriver.entrySet()) {
            EnrichedSomaticVariant variant = entry.getKey();
            if (!variantsToReport.contains(variant) && !Collections.disjoint(entry.getValue(), filteredEvidence)) {
                LOGGER.warn("Evidence found on somatic variant on gene {} which is not included in driver catalog!", variant.gene());
            }
        }

        // Check that we miss no drivers
        for (DriverCatalog driver : driverCatalog) {
            boolean reported = false;
            for (EnrichedSomaticVariant variant : variantsToReport) {
                if (variant.gene().equals(driver.gene())) {
                    reported = true;
                }
            }
            if (!reported) {
                LOGGER.warn("Driver {} not added to reported somatic variants", driver.gene());
            }
        }

        return ImmutableSomaticVariantAnalysis.of(variantsToReport,
                filteredEvidence,
                driverCatalog,
                microsatelliteIndelsPerMb,
                tumorMutationalLoad,
                tumorMutationalBurden);
    }

    @NotNull
    private static Predicate<EnrichedSomaticVariant> includeFilter(@NotNull GeneModel panelGeneModel) {
        return variant -> {
            if (variant.isFiltered()) {
                return false;
            }

            if (!panelGeneModel.somaticVariantGenes().contains(variant.gene())) {
                return false;
            }

            CodingEffect effect = variant.canonicalCodingEffect();
            DriverCategory category = panelGeneModel.category(variant.gene());

            if (category == null) {
                throw new IllegalStateException("Driver category not known for driver gene: " + variant.gene());
            }

            if (variant.isHotspot()) {
                return true;
            } else if (category == DriverCategory.TSG) {
                return TSG_CODING_EFFECTS_TO_REPORT.contains(effect);
            } else {
                assert category == DriverCategory.ONCO;
                if (ONCO_GENES_WITH_SPLICE_EFFECTS.contains(variant.gene())) {
                    return ONCO_CODING_EFFECTS_TO_REPORT.contains(effect) || effect == CodingEffect.SPLICE;
                } else {
                    return ONCO_CODING_EFFECTS_TO_REPORT.contains(effect);
                }
            }
        };
    }
}
