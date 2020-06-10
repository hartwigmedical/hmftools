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
import com.hartwig.hmftools.serve.vicc.range.GeneLevelEventExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeExtractor;
import com.hartwig.hmftools.serve.vicc.signatures.SignaturesExtractor;
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
    @NotNull
    private final GeneLevelEventExtractor geneLevelEventExtractor;
    @NotNull
    private final GeneRangeExtractor geneRangeExtractor;
    @NotNull
    private final SignaturesExtractor signaturesExtractor;

    public ViccExtractor(@NotNull final HotspotExtractor hotspotExtractor, @NotNull final CopyNumberExtractor copyNumberExtractor,
            @NotNull final FusionExtractor fusionExtractor, @NotNull final GeneLevelEventExtractor geneLevelEventExtractor,
            @NotNull final GeneRangeExtractor geneRangeExtractor, @NotNull final SignaturesExtractor signaturesExtractor) {
        this.hotspotExtractor = hotspotExtractor;
        this.copyNumberExtractor = copyNumberExtractor;
        this.fusionExtractor = fusionExtractor;
        this.geneLevelEventExtractor = geneLevelEventExtractor;
        this.geneRangeExtractor = geneRangeExtractor;
        this.signaturesExtractor = signaturesExtractor;
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
            Map<Feature, String> geneLevelEventsPerFeature = geneLevelEventExtractor.extractKnownGeneLevelEvents(viccEntry);
            Map<Feature, String> geneRangesPerFeature = geneRangeExtractor.extractGeneRanges(viccEntry);
            Map<Feature, String> signaturesPerFeature = signaturesExtractor.extractSignatures(viccEntry);

            extractionResultsPerEntry.put(viccEntry,
                    ImmutableViccExtractionResult.builder()
                            .hotspotsPerFeature(hotspotsPerFeature)
                            .ampsDelsPerFeature(ampsDelsPerFeature)
                            .fusionsPerFeature(fusionsPerFeature)
                            .geneLevelEventsPerFeature(geneLevelEventsPerFeature)
                            .geneRangesPerFeature(geneRangesPerFeature)
                            .signaturesPerFeature(signaturesPerFeature)
                            .build());
        }

        LOGGER.info("Could not resolve hotspots for {} features", hotspotExtractor.unresolvableFeatures().size());
        for (String feature : hotspotExtractor.unresolvableFeatures()) {
           // LOGGER.debug(" {}", feature);
        }

        return extractionResultsPerEntry;
    }
}
