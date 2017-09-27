package com.hartwig.hmftools.patientreporter.report.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.apiclients.civic.api.CivicApiWrapper;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;
import com.hartwig.hmftools.apiclients.diseaseontology.api.DiseaseOntologyApiWrapper;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.jasperreports.engine.data.JRBeanCollectionDataSource;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class EvidenceReportData {
    private static final Logger LOGGER = LogManager.getLogger(EvidenceReportData.class);

    @NotNull
    public abstract List<Alteration> getAlterations();

    @NotNull
    @Value.Derived
    public List<Alteration> getAlterationsWithEvidence() {
        return getAlterations().stream().filter(alteration -> alteration.getEvidence().size() > 0).collect(Collectors.toList());
    }

    public static EvidenceReportData of(@NotNull final PatientReport report, @NotNull final GeneModel geneModel,
            @NotNull final Set<String> tumorDoids) {
        final List<Alteration> alterations = Lists.newArrayList();
        final Set<String> tumorSubtypesDoids = getTumorSubtypesDoids(tumorDoids);
        for (final VariantReport variantReport : report.variants()) {
            for (final HmfGenomeRegion region : geneModel.hmfRegions()) {
                if (region.gene().equals(variantReport.gene())) {
                    final int entrezId = Integer.parseInt(region.entrezId());
                    final List<CivicVariant> civicVariants = getCivicVariants(entrezId, variantReport);
                    final Alteration alteration = Alteration.from(variantReport, civicVariants, tumorSubtypesDoids);
                    if (alteration.getMatches().size() > 0) {
                        alterations.add(alteration);
                    }
                }
            }
        }
        CivicApiWrapper.releaseResources();
        DiseaseOntologyApiWrapper.releaseResources();
        return ImmutableEvidenceReportData.of(alterations);
    }

    private static Set<String> getTumorSubtypesDoids(@NotNull final Set<String> tumorDoids) {
        final Set<String> tumorSubtypesDoids = Sets.newHashSet();
        tumorSubtypesDoids.addAll(tumorDoids);
        for (final String tumorDoid : tumorDoids) {
            try {
                final List<String> childrenDoid = DiseaseOntologyApiWrapper.getAllChildrenDoids(tumorDoid).toList().blockingGet();
                tumorSubtypesDoids.addAll(childrenDoid);
            } catch (final Throwable throwable) {
                LOGGER.error("Failed to get children doids for tumor doid: " + tumorDoid + ". error message: " + throwable.getMessage());
            }
        }
        return tumorSubtypesDoids;
    }

    @NotNull
    private static List<CivicVariant> getCivicVariants(final int entrezId, @NotNull final VariantReport variantReport) {
        try {
            return CivicApiWrapper.getVariantsForGene(entrezId).toList().blockingGet();
        } catch (final Throwable throwable) {
            LOGGER.error("Failed to get civic variants for variant: " + variantReport.variant().chromosomePosition() + ". error message: "
                    + throwable.getMessage());
            return Lists.newArrayList();
        }
    }

    @NotNull
    public JRBeanCollectionDataSource toDataSource() {
        return new JRBeanCollectionDataSource(Lists.newArrayList(this));
    }
}
