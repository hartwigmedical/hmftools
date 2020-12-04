package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class HotspotExtractor {

    @NotNull
    private final GeneChecker geneChecker;
    @NotNull
    private final ProteinResolver proteinResolver;
    @NotNull
    private final ProteinAnnotationExtractor proteinAnnotationExtractor;

    public HotspotExtractor(@NotNull final GeneChecker geneChecker, @NotNull final ProteinResolver proteinResolver,
            @NotNull final ProteinAnnotationExtractor proteinAnnotationExtractor) {
        this.geneChecker = geneChecker;
        this.proteinResolver = proteinResolver;
        this.proteinAnnotationExtractor = proteinAnnotationExtractor;
    }

    @NotNull
    public Map<Feature, List<VariantHotspot>> extractHotspots(@NotNull ViccEntry viccEntry) {
        Map<Feature, List<VariantHotspot>> hotspotsPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            if (feature.type() == MutationType.HOTSPOT && geneChecker.isValidGene(feature.geneSymbol(), feature.name())) {
                String proteinAnnotation = proteinAnnotationExtractor.apply(feature.name());
                List<VariantHotspot> hotspots = proteinResolver.resolve(feature.geneSymbol(), viccEntry.transcriptId(), proteinAnnotation);
                hotspotsPerFeature.put(feature, hotspots);
            }
        }

        return hotspotsPerFeature;
    }

    @NotNull
    public ProteinAnnotationExtractor proteinAnnotationExtractor() {
        return proteinAnnotationExtractor;
    }
}
