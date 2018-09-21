package com.hartwig.hmftools.common.variant.snpeff;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.region.CanonicalTranscriptFactory;
import com.hartwig.hmftools.common.region.TranscriptRegion;
import com.hartwig.hmftools.common.variant.cosmic.CosmicAnnotation;

import org.jetbrains.annotations.NotNull;

public class CanonicalAnnotation {

    @NotNull
    private final Set<String> driverCatalogGenes;
    @NotNull
    private final Map<String, String> canonicalTranscriptGeneMap;

    public CanonicalAnnotation() {
        this.driverCatalogGenes = asSet(DndsDriverGeneLikelihoodSupplier.tsgLikelihood());
        this.driverCatalogGenes.addAll(asSet(DndsDriverGeneLikelihoodSupplier.oncoLikelihood()));
        this.canonicalTranscriptGeneMap = CanonicalTranscriptFactory.create()
                .stream()
                .collect(Collectors.toMap(TranscriptRegion::transcriptID, TranscriptRegion::gene));
    }

    @NotNull
    public Optional<CosmicAnnotation> canonicalCosmicAnnotation(@NotNull final List<CosmicAnnotation> cosmicAnnotations) {
        return cosmicAnnotations.stream().filter(x -> canonicalTranscriptGeneMap.keySet().contains(x.transcript())).findFirst();
    }

    @NotNull
    public Optional<SnpEffAnnotation> canonicalSnpEffAnnotation(@NotNull final List<SnpEffAnnotation> allAnnotations) {
        final List<SnpEffAnnotation> canonicalTranscriptAnnotations = allAnnotations.stream()
                .filter(SnpEffAnnotation::isTranscriptFeature)
                .filter(x -> canonicalTranscriptGeneMap.keySet().contains(x.transcript()))
                .collect(Collectors.toList());

        if (!canonicalTranscriptAnnotations.isEmpty()) {
            final Optional<SnpEffAnnotation> canonicalOnDriverGene =
                    canonicalTranscriptAnnotations.stream().filter(x -> driverCatalogGenes.contains(x.gene())).findFirst();
            if (canonicalOnDriverGene.isPresent()) {
                return canonicalOnDriverGene;
            }

            return Optional.of(canonicalTranscriptAnnotations.get(0));
        }

        return Optional.empty();
    }

    @NotNull
    private static Set<String> asSet(Map<String, DndsDriverGeneLikelihood> dndsLikelihoods) {
        return dndsLikelihoods.values().stream().map(DndsDriverGeneLikelihood::gene).collect(Collectors.toSet());
    }
}
