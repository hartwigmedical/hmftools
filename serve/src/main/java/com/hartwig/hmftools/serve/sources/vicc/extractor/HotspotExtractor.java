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

import org.jetbrains.annotations.NotNull;

public class HotspotExtractor {

    @NotNull
    private final ProteinResolver proteinResolver;

    public HotspotExtractor(@NotNull final ProteinResolver proteinResolver) {
        this.proteinResolver = proteinResolver;
    }

    @NotNull
    public Map<Feature, List<VariantHotspot>> extractHotspots(@NotNull ViccEntry viccEntry) {
        Map<Feature, List<VariantHotspot>> hotspotsPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            if (feature.type() == MutationType.HOTSPOT) {
                hotspotsPerFeature.put(feature,
                        proteinResolver.extractHotspotsFromProteinAnnotation(feature.geneSymbol(),
                                viccEntry.transcriptId(),
                                ProteinAnnotationExtractor.proteinAnnotation(feature)));
            }
        }

        return hotspotsPerFeature;
    }
}
