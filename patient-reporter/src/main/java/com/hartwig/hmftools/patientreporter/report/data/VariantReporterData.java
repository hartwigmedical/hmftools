package com.hartwig.hmftools.patientreporter.report.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.apiclients.civic.api.CivicApiWrapper;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;
import com.hartwig.hmftools.apiclients.civic.data.ImmutableCivicVariant;
import com.hartwig.hmftools.apiclients.diseaseontology.api.DiseaseOntologyApiWrapper;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.jasperreports.engine.data.JRBeanCollectionDataSource;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class VariantReporterData {
    private static final Logger LOGGER = LogManager.getLogger(VariantReporterData.class);

    @NotNull
    public abstract List<AlterationReporterData> getAlterations();

    public static VariantReporterData of(@NotNull final PatientReport report, @NotNull final GeneModel geneModel,
            @NotNull final Set<String> tumorDoids) {
        final List<AlterationReporterData> alterations = Lists.newArrayList();
        final Set<String> tumorChildrenDoids = getTumorChildrenDoids(tumorDoids);
        for (final VariantReport variantReport : report.variants()) {
            for (final HmfGenomeRegion region : geneModel.hmfRegions()) {
                if (region.gene().equals(variantReport.gene())) {
                    final int entrezId = Integer.parseInt(region.entrezId());
                    final List<CivicVariant> civicVariants = getCivicVariants(entrezId, variantReport, tumorChildrenDoids);
                    final AlterationReporterData alteration =
                            AlterationReporterData.from(variantReport.gene(), variantReport.hgvsProtein(), "", "", civicVariants);
                    if (alteration.getEvidence().size() > 0) {
                        alterations.add(alteration);
                    }
                }
            }
        }
        CivicApiWrapper.releaseResources();
        DiseaseOntologyApiWrapper.releaseResources();
        return ImmutableVariantReporterData.of(alterations);
    }

    private static Set<String> getTumorChildrenDoids(@NotNull final Set<String> tumorDoids) {
        final Set<String> tumorChildrenDoids = Sets.newHashSet();
        tumorChildrenDoids.addAll(tumorDoids);
        for (final String tumorDoid : tumorDoids) {
            try {
                final List<String> childrenDoid = DiseaseOntologyApiWrapper.getAllChildrenDoids(tumorDoid).toList().blockingGet();
                tumorChildrenDoids.addAll(childrenDoid);
            } catch (final Throwable throwable) {
                LOGGER.error("Failed to get children doids for tumor doid: " + tumorDoid + ". error message: " + throwable.getMessage());
            }
        }
        return tumorChildrenDoids;
    }

    @NotNull
    private static List<CivicVariant> getCivicVariants(final int entrezId, @NotNull final VariantReport variantReport,
            @NotNull final Set<String> tumorChildrenDoids) {
        final Variant variant = variantReportToVariant(variantReport);
        try {
            //MIVO: get exact matches only for now (will not work for indels)
            return CivicApiWrapper.getVariantMatches(entrezId, variant)
                    .map(civicVariant -> (CivicVariant) ImmutableCivicVariant.builder()
                            .from(civicVariant)
                            .evidenceItems(civicVariant.evidenceItemsWithDrugs()
                                    .stream()
                                    .filter(evidenceItem -> tumorChildrenDoids.contains(evidenceItem.disease().doidString()))
                                    .collect(Collectors.toList()))
                            .build())
                    .filter(civicVariant -> !civicVariant.groupedEvidenceItems().isEmpty())
                    .toList()
                    .blockingGet();
        } catch (final Throwable throwable) {
            LOGGER.error("Failed to get civic variants for variant: " + variant.chromosomePosition() + ". error message: "
                    + throwable.getMessage());
            return Lists.newArrayList();
        }
    }

    @NotNull
    public JRBeanCollectionDataSource toDataSource() {
        if (getAlterations().size() > 0) {
            return new JRBeanCollectionDataSource(Lists.newArrayList(this));
        } else {
            return new JRBeanCollectionDataSource(Lists.newArrayList());
        }
    }

    public static Variant variantReportToVariant(@NotNull final VariantReport variantReport) {
        return new Variant() {
            @NotNull
            @Override
            public String ref() {
                return variantReport.ref();
            }

            @NotNull
            @Override
            public String alt() {
                return variantReport.alt();
            }

            @NotNull
            @Override
            public VariantType type() {
                return null;
            }

            @NotNull
            @Override
            public String filter() {
                return null;
            }

            @NotNull
            @Override
            public String chromosome() {
                return variantReport.chromosome();
            }

            @Override
            public long position() {
                return variantReport.position();
            }
        };
    }
}
