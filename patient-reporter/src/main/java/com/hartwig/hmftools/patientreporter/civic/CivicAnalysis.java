package com.hartwig.hmftools.patientreporter.civic;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.apiclients.civic.api.CivicApiWrapper;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;
import com.hartwig.hmftools.apiclients.diseaseontology.api.DiseaseOntologyApiWrapper;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.patientreporter.report.data.Alteration;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CivicAnalysis {
    private static final Logger LOGGER = LogManager.getLogger(CivicAnalysis.class);

    @NotNull
    public static List<Alteration> run(@NotNull final List<VariantReport> reportedVariants, @NotNull final List<GeneCopyNumber> copyNumbers,
            @NotNull final GeneModel geneModel, @NotNull final Set<String> tumorDoids) {
        LOGGER.info(" Analysing civic associations...");
        final Set<String> relevantDoids = getRelevantDoids(tumorDoids);
        if (relevantDoids.isEmpty()) {
            LOGGER.warn("  Disease-ontology id set for this tumor is empty!");
        }
        final List<Alteration> alterations = civicVariantAlterations(reportedVariants, geneModel, relevantDoids);
        alterations.addAll(civicCopyNumberAlterations(copyNumbers, geneModel, relevantDoids));
        return alterations;
    }

    @NotNull
    private static Set<String> getRelevantDoids(@NotNull final Set<String> tumorDoids) {
        LOGGER.info("  Fetching relevant tumor types...");
        final Set<String> relevantDoids = Sets.newHashSet();
        relevantDoids.addAll(tumorDoids);
        final DiseaseOntologyApiWrapper diseaseOntologyApi = new DiseaseOntologyApiWrapper();
        for (final String tumorDoid : tumorDoids) {
            try {
                final List<String> childrenAndParentDoids = diseaseOntologyApi.getAllChildrenDoids(tumorDoid)
                        .mergeWith(diseaseOntologyApi.getAllParentDoids(tumorDoid))
                        .toList()
                        .blockingGet();
                relevantDoids.addAll(childrenAndParentDoids);
            } catch (final Throwable throwable) {
                LOGGER.error("  Failed to get children doids for tumor doid: " + tumorDoid + ". error message: " + throwable.getMessage());
            }
        }
        diseaseOntologyApi.releaseResources();
        return relevantDoids;
    }

    @NotNull
    private static List<Alteration> civicVariantAlterations(@NotNull final List<VariantReport> reportedVariants,
            @NotNull final GeneModel geneModel, @NotNull final Set<String> relevantDoids) {
        LOGGER.info("  Fetching civic variant alterations...");
        final List<Alteration> alterations = Lists.newArrayList();
        final CivicApiWrapper civicApi = new CivicApiWrapper();
        for (final VariantReport variantReport : reportedVariants) {
            for (final HmfGenomeRegion region : geneModel.regions()) {
                if (region.gene().equals(variantReport.gene())) {
                    for (final int entrezId : region.entrezId()) {
                        try {
                            final List<CivicVariant> civicVariants = civicApi.getVariantsForGene(entrezId).toList().blockingGet();
                            final Alteration alteration = Alteration.from(variantReport, civicVariants, relevantDoids);
                            if (alteration.getMatches().size() > 0) {
                                alterations.add(alteration);
                            }
                        } catch (final Throwable throwable) {
                            LOGGER.error("  Failed to get civic variants for variant: " + variantReport.variant().chromosomePosition()
                                    + ". error message: " + throwable.getMessage());
                        }
                    }
                }
            }
        }
        civicApi.releaseResources();
        return alterations;
    }

    @NotNull
    private static List<Alteration> civicCopyNumberAlterations(@NotNull final List<GeneCopyNumber> copyNumbers,
            @NotNull final GeneModel geneModel, @NotNull final Set<String> relevantDoids) {
        LOGGER.info("  Fetching civic copy number alterations...");
        final List<Alteration> alterations = Lists.newArrayList();
        final CivicApiWrapper civicApi = new CivicApiWrapper();
        for (final GeneCopyNumber copyNumberReport : copyNumbers) {
            for (final HmfGenomeRegion region : geneModel.regions()) {
                if (region.gene().equals(copyNumberReport.gene())) {
                    for (final int entrezId : region.entrezId()) {
                        try {
                            final List<CivicVariant> civicVariants = civicApi.getVariantsForGene(entrezId).toList().blockingGet();
                            final Alteration alteration = Alteration.from(copyNumberReport, civicVariants, relevantDoids);
                            if (alteration.getMatches().size() > 0) {
                                alterations.add(alteration);
                            }
                        } catch (final Throwable throwable) {
                            LOGGER.error("  Failed to get civic variants for copy number: {}. error message: {}",
                                    copyNumberReport.gene(),
                                    throwable.getMessage());
                        }
                    }
                }
            }
        }
        civicApi.releaseResources();
        return alterations;
    }
}
