package com.hartwig.hmftools.serve.vicc;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberAnnotation;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.vicc.fusion.FusionAnnotation;
import com.hartwig.hmftools.serve.vicc.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.vicc.genelevel.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.vicc.genelevel.GeneLevelEventExtractor;
import com.hartwig.hmftools.serve.vicc.hotspot.HotspotExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeAnnotation;
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
    public Map<ViccEntry, ViccExtractionResult> extractFromViccEntries(@NotNull List<ViccEntry> viccEntries) {
        Map<ViccEntry, ViccExtractionResult> extractionResultsPerEntry = Maps.newHashMap();
        for (ViccEntry entry : viccEntries) {
            Map<Feature, List<VariantHotspot>> hotspotsPerFeature = hotspotExtractor.extractHotspots(entry);
            Map<Feature, CopyNumberAnnotation> ampsDelsPerFeature = copyNumberExtractor.extractKnownAmplificationsDeletions(entry);
            Map<Feature, FusionAnnotation> fusionsPerFeature = fusionExtractor.extractKnownFusions(entry);
            Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = geneLevelEventExtractor.extractKnownGeneLevelEvents(entry);
            Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = geneRangeExtractor.extractGeneRanges(entry);
            Map<Feature, SignatureName> signaturesPerFeature = signaturesExtractor.extractSignatures(entry);

            ActionableEvidence actionableEvidence = ActionableEvidenceFactory.toActionableEvidence(entry);

            extractionResultsPerEntry.put(entry,
                    ImmutableViccExtractionResult.builder()
                            .hotspotsPerFeature(hotspotsPerFeature)
                            .ampsDelsPerFeature(ampsDelsPerFeature)
                            .fusionsPerFeature(fusionsPerFeature)
                            .geneLevelEventsPerFeature(geneLevelEventsPerFeature)
                            .geneRangesPerFeature(geneRangesPerFeature)
                            .signaturesPerFeature(signaturesPerFeature)
                            .actionableEvidence(actionableEvidence)
                            .build());
        }

        LOGGER.info("Unique known amps: {}", copyNumberExtractor.uniqueAmps().size());
        LOGGER.info("Unique known dels: {}", copyNumberExtractor.uniqueDels().size());
        LOGGER.info("Unique known fusion pairs: {}", fusionExtractor.uniqueFusionsPair().size());

        return extractionResultsPerEntry;
    }


}
