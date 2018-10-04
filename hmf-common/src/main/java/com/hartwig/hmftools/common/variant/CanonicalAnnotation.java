package com.hartwig.hmftools.common.variant;

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
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import org.jetbrains.annotations.NotNull;

class CanonicalAnnotation {

    @NotNull
    private final Set<String> driverCatalogGenes;
    @NotNull
    private final Map<String, String> canonicalTranscriptGeneMap;

    CanonicalAnnotation() {
        this.driverCatalogGenes = asSet(DndsDriverGeneLikelihoodSupplier.tsgLikelihood());
        this.driverCatalogGenes.addAll(asSet(DndsDriverGeneLikelihoodSupplier.oncoLikelihood()));
        this.canonicalTranscriptGeneMap = CanonicalTranscriptFactory.create()
                .stream()
                .collect(Collectors.toMap(TranscriptRegion::transcriptID, TranscriptRegion::gene));
    }

    @NotNull
    public Optional<CosmicAnnotation> canonicalCosmicAnnotation(@NotNull final List<CosmicAnnotation> cosmicAnnotations) {
        return pickCanonicalFavourDriverGene(cosmicAnnotations);
    }

    @NotNull
    public Optional<SnpEffAnnotation> canonicalSnpEffAnnotation(@NotNull final List<SnpEffAnnotation> allAnnotations) {
        final List<SnpEffAnnotation> transcriptAnnotations = allAnnotations.stream()
                .filter(SnpEffAnnotation::isTranscriptFeature).collect(Collectors.toList());
        return pickCanonicalFavourDriverGene(transcriptAnnotations);
    }

    @NotNull
    <T extends TranscriptAnnotation> Optional<T> pickCanonicalFavourDriverGene(@NotNull List<T> annotations) {
        final List<T> canonicalAnnotations = annotations.stream()
                .filter(annotation -> canonicalTranscriptGeneMap.keySet().contains(annotation.transcript()))
                .collect(Collectors.toList());

        if (!canonicalAnnotations.isEmpty()) {
            final Optional<T> canonicalOnDriverGene =
                    canonicalAnnotations.stream().filter(annotation -> driverCatalogGenes.contains(annotation.gene())).findFirst();
            if (canonicalOnDriverGene.isPresent()) {
                return canonicalOnDriverGene;
            }

            return Optional.of(canonicalAnnotations.get(0));
        }

        return Optional.empty();
    }

    @NotNull
    private static Set<String> asSet(Map<String, DndsDriverGeneLikelihood> dndsLikelihoods) {
        return dndsLikelihoods.values().stream().map(DndsDriverGeneLikelihood::gene).collect(Collectors.toSet());
    }
}
