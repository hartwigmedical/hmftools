package com.hartwig.hmftools.serve.vicc.hotspot;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.HotspotGenerator;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class HotspotExtractor {

    @NotNull
    private final HotspotGenerator hotspotGenerator;

    public HotspotExtractor(@NotNull final HotspotGenerator hotspotGenerator) {
        this.hotspotGenerator = hotspotGenerator;
    }

    @NotNull
    public Map<Feature, List<VariantHotspot>> extractHotspots(@NotNull ViccEntry viccEntry) {
        Map<Feature, List<VariantHotspot>> allHotspotsPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            String proteinAnnotation = feature.proteinAnnotation();
            if (HotspotGenerator.isResolvableProteinAnnotation(proteinAnnotation)) {
                if (proteinAnnotation.equals("S257delinsK") || proteinAnnotation.equals("*214W")) {
                    int x = 1;
                }
                allHotspotsPerFeature.put(feature,
                        hotspotGenerator.generateHotspots(feature.geneSymbol(), viccEntry.transcriptId(), proteinAnnotation));
            }
        }

        return allHotspotsPerFeature;
    }

}
