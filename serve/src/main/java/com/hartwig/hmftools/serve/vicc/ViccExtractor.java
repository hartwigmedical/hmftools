package com.hartwig.hmftools.serve.vicc;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.vicc.hotspot.HotspotExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ViccExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ViccExtractor.class);

    @NotNull
    private final HotspotExtractor hotspotExtractor;
    @NotNull
    private final CopyNumberExtractor copyNumberExtractor;
    @NotNull
    private final FusionExtractor fusionExtractor;

    public ViccExtractor(@NotNull final HotspotExtractor hotspotExtractor, @NotNull final CopyNumberExtractor copyNumberExtractor,
            @NotNull final FusionExtractor fusionExtractor) {
        this.hotspotExtractor = hotspotExtractor;
        this.copyNumberExtractor = copyNumberExtractor;
        this.fusionExtractor = fusionExtractor;
    }

    @NotNull
    public Map<ViccEntry, ViccExtractionResult> extractFromViccEntries(@NotNull List<ViccEntry> viccEntries)
            throws IOException, InterruptedException {
        Map<ViccEntry, ViccExtractionResult> extractionResultsPerEntry = Maps.newHashMap();
        for (ViccEntry viccEntry : viccEntries) {
            Map<Feature, List<VariantHotspot>> hotspotsPerFeature = hotspotExtractor.extractHotspots(viccEntry);
            Map<Feature, KnownAmplificationDeletion> ampsDelsPerFeature =
                    copyNumberExtractor.extractKnownAmplificationsDeletions(viccEntry);
            Map<Feature, String> fusionsPerFeature = fusionExtractor.extractKnownFusions(viccEntry);

            extractionResultsPerEntry.put(viccEntry,
                    ImmutableViccExtractionResult.builder()
                            .hotspotsPerFeature(hotspotsPerFeature)
                            .ampsDelsPerFeature(ampsDelsPerFeature)
                            .fusionsPerFeature(fusionsPerFeature)
                            .build());
        }

        LOGGER.info("Could not resolve hotspots for {} features", hotspotExtractor.unresolvableFeatures().size());
        for (String feature : hotspotExtractor.unresolvableFeatures()) {
            LOGGER.debug(" {}", feature);
        }

        return extractionResultsPerEntry;
    }
}
