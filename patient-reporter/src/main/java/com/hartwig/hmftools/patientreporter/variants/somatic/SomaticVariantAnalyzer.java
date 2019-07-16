package com.hartwig.hmftools.patientreporter.variants.somatic;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
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
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.variants.driver.DriverGeneView;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SomaticVariantAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantAnalyzer.class);

    private static final List<CodingEffect> TSG_CODING_EFFECTS_TO_REPORT =
            Lists.newArrayList(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE, CodingEffect.SPLICE);

    private static final Set<String> ONCO_GENES_WITH_SPLICE_EFFECTS = Sets.newHashSet("MET", "JAK2");

    private static final List<CodingEffect> ONCO_CODING_EFFECTS_TO_REPORT = Lists.newArrayList(CodingEffect.MISSENSE);

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static SomaticVariantAnalysis run(@NotNull List<EnrichedSomaticVariant> variants, @NotNull DriverGeneView driverGeneView,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation) {
        List<EnrichedSomaticVariant> variantsToReport =
                variants.stream().filter(includeFilter(driverGeneView)).collect(Collectors.toList());

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant =
                actionabilityAnalyzer.evidenceForSomaticVariants(variants, primaryTumorLocation);

        List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(OncoDrivers.drivers(variants));
        driverCatalog.addAll(TsgDrivers.drivers(variants));

        // Extract somatic evidence for high drivers variants into flat list (See DEV-824)
        List<EvidenceItem> filteredEvidence =
                ReportableEvidenceItemFactory.reportableFlatListDriversOnly(evidencePerVariant, driverCatalog);

        // Check that all variants with high level evidence are reported (since they are in the driver catalog).
        for (Map.Entry<EnrichedSomaticVariant, List<EvidenceItem>> entry : evidencePerVariant.entrySet()) {
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
                MicrosatelliteAnalyzer.determineMicrosatelliteIndelsPerMb(variants),
                MutationalLoadAnalyzer.determineTumorMutationalLoad(variants),
                MutationalBurdenAnalyzer.determineTumorMutationalBurden(variants));
    }

    @NotNull
    private static Predicate<EnrichedSomaticVariant> includeFilter(@NotNull DriverGeneView driverGeneView) {
        return variant -> {
            if (variant.isFiltered()) {
                return false;
            }

            if (!driverGeneView.oncoDriverGenes().contains(variant.gene()) && !driverGeneView.tsgDriverGenes().contains(variant.gene())) {
                return false;
            }

            // Report all hotspots on our driver genes.
            if (variant.isHotspot()) {
                return true;
            }

            DriverCategory category = driverGeneView.category(variant.gene());
            if (category == null) {
                throw new IllegalStateException("Driver category not known for driver gene: " + variant.gene());
            }

            CodingEffect effect = variant.canonicalCodingEffect();
            if (category == DriverCategory.TSG) {
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
