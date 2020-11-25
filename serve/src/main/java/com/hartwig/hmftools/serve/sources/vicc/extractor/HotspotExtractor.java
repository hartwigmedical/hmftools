package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class HotspotExtractor {

    @NotNull
    private final ProteinResolver proteinResolver;
    @NotNull
    private final ProteinAnnotationExtractor proteinAnnotationExtractor;

    public HotspotExtractor(@NotNull final ProteinResolver proteinResolver,
            @NotNull final ProteinAnnotationExtractor proteinAnnotationExtractor) {
        this.proteinResolver = proteinResolver;
        this.proteinAnnotationExtractor = proteinAnnotationExtractor;
    }

    @NotNull
    public Map<Feature, List<VariantHotspot>> extractHotspots(@NotNull ViccEntry viccEntry) {
        Map<Feature, List<VariantHotspot>> hotspotsPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            if (feature.type() == MutationType.HOTSPOT) {
                hotspotsPerFeature.put(feature,
                        proteinResolver.extractHotspotsFromProteinAnnotation(feature.geneSymbol(),
                                viccEntry.transcriptId(),
                                extractProteinAnnotation(feature)));
            }
        }

        return hotspotsPerFeature;
    }

    @NotNull
    public String extractProteinAnnotation(@NotNull Feature feature) {
        return feature.type() == MutationType.HOTSPOT ? proteinAnnotationExtractor.apply(feature.name()) : Strings.EMPTY;
    }
}
