package com.hartwig.hmftools.serve.vicc.hotspot;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinToHotspotConverter;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class HotspotExtractor {

    @NotNull
    private final ProteinToHotspotConverter proteinToHotspotConverter;

    public HotspotExtractor(@NotNull final ProteinToHotspotConverter proteinToHotspotConverter) {
        this.proteinToHotspotConverter = proteinToHotspotConverter;
    }

    @NotNull
    public Map<Feature, List<VariantHotspot>> extractHotspots(@NotNull ViccEntry viccEntry) {
        Map<Feature, List<VariantHotspot>> allHotspotsPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            String proteinAnnotation = feature.proteinAnnotation();
            if (ProteinToHotspotConverter.isResolvableProteinAnnotation(proteinAnnotation)) {
                allHotspotsPerFeature.put(feature,
                        proteinToHotspotConverter.resolveProteinAnnotation(feature.geneSymbol(), viccEntry.transcriptId(), proteinAnnotation));
            }
        }

        return allHotspotsPerFeature;
    }

}
