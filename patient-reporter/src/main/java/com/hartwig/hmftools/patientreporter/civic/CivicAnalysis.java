package com.hartwig.hmftools.patientreporter.civic;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Predicate;
import java.util.stream.Collectors;

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

import io.reactivex.Observable;

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
        final Collection<HmfGenomeRegion> geneRegions = geneModel.regions().stream().filter(distinctByGene()).collect(Collectors.toList());
        final List<Alteration> alterations = civicVariantAlterations(reportedVariants, geneRegions, relevantDoids);
        alterations.addAll(civicCopyNumberAlterations(copyNumbers, geneRegions, relevantDoids));
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
            @NotNull final Collection<HmfGenomeRegion> geneRegions, @NotNull final Set<String> relevantDoids) {
        LOGGER.info("  Fetching civic variant alterations...");
        final List<Alteration> alterations = Lists.newArrayList();
        final CivicApiWrapper civicApi = new CivicApiWrapper();
        final List<HmfGenomeRegion> wildTypeGenes = Lists.newArrayList();
        for (final HmfGenomeRegion region : geneRegions) {
            final List<VariantReport> reportedVariantsInGene = reportedVariants.stream()
                    .filter(variantReport -> region.gene().equals(variantReport.gene()))
                    .collect(Collectors.toList());
            if (!reportedVariantsInGene.isEmpty()) {
                for (final VariantReport variantReport : reportedVariantsInGene) {
                    alterations.addAll(queryCivicAlteration(region,
                            variantList -> Alteration.from(variantReport, variantList, relevantDoids),
                            "  Failed to get civic variants for variant: " + variantReport.variant().chromosomePosition()));
                }
            } else {
                wildTypeGenes.add(region);
            }
        }
        final List<Alteration> wildTypeAlterations = Observable.fromIterable(wildTypeGenes)
                .flatMap(hmfGenomeRegion -> civicApi.getWildTypeVariantsForGeneIds(hmfGenomeRegion.entrezId())
                        .toList()
                        .map(variantList -> Alteration.fromWildType(hmfGenomeRegion.gene(), variantList, relevantDoids))
                        .filter(alteration -> !alteration.getMatches().isEmpty())
                        .toObservable())
                .toList()
                .blockingGet();
        alterations.addAll(wildTypeAlterations);
        civicApi.releaseResources();
        return alterations;
    }

    @NotNull
    private static List<Alteration> civicCopyNumberAlterations(@NotNull final List<GeneCopyNumber> copyNumbers,
            @NotNull final Collection<HmfGenomeRegion> geneRegions, @NotNull final Set<String> relevantDoids) {
        LOGGER.info("  Fetching civic copy number alterations...");
        final List<Alteration> alterations = Lists.newArrayList();
        for (final GeneCopyNumber copyNumberReport : copyNumbers) {
            for (final HmfGenomeRegion region : geneRegions) {
                if (region.gene().equals(copyNumberReport.gene())) {
                    alterations.addAll(queryCivicAlteration(region,
                            variantList -> Alteration.from(copyNumberReport, variantList, relevantDoids),
                            "  Failed to get civic variants for copy number: " + copyNumberReport.gene()));
                }
            }
        }
        return alterations;
    }

    @NotNull
    private static List<Alteration> queryCivicAlteration(@NotNull final HmfGenomeRegion region,
            io.reactivex.functions.Function<List<CivicVariant>, Alteration> alterationBuilder, @NotNull final String error) {
        final CivicApiWrapper civicApi = new CivicApiWrapper();
        final List<Alteration> result = Lists.newArrayList();
        try {
            result.addAll(Observable.fromIterable(region.entrezId())
                    .flatMap(civicApi::getVariantsForGene)
                    .toList()
                    .map(alterationBuilder)
                    .filter(alteration -> !alteration.getMatches().isEmpty())
                    .toObservable()
                    .toList()
                    .blockingGet());
        } catch (final Throwable throwable) {
            LOGGER.error("{}. Error message: {}", error, throwable.getMessage());
        }
        civicApi.releaseResources();
        return result;
    }

    @NotNull
    private static Predicate<HmfGenomeRegion> distinctByGene() {
        final Set<String> seen = ConcurrentHashMap.newKeySet();
        return t -> seen.add(t.geneID());
    }
}
