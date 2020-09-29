package com.hartwig.hmftools.serve.vicc;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.fusion.FusionAnnotation;
import com.hartwig.hmftools.serve.vicc.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.vicc.hotspot.HotspotExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneLevelEventExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeAnnotation;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeExtractor;
import com.hartwig.hmftools.serve.vicc.signatures.SignaturesExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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
            Map<Feature, KnownAmplificationDeletion> ampsDelsPerFeature = copyNumberExtractor.extractKnownAmplificationsDeletions(entry);
            Map<Feature, FusionAnnotation> fusionsPerFeature = fusionExtractor.extractKnownFusions(entry);
            Map<Feature, String> geneLevelEventsPerFeature = geneLevelEventExtractor.extractKnownGeneLevelEvents(entry);
            Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = geneRangeExtractor.extractGeneRanges(entry);
            Map<Feature, String> signaturesPerFeature = signaturesExtractor.extractSignatures(entry);

            ActionableEvidence actionableEvidence = extractEvidence(entry);

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
        LOGGER.info("Unique known fusion promiscuous: {}", fusionExtractor.uniqueFusionsPromiscuous().size());

        return extractionResultsPerEntry;
    }

    @Nullable
    private static ActionableEvidence extractEvidence(@NotNull ViccEntry entry) {
        String drugs = entry.association().drugLabels();

        String cancerTypeString = null;
        String cancerTypeDOID = null;
        Phenotype phenotype = entry.association().phenotype();
        if (phenotype != null) {
            cancerTypeString = phenotype.description();
            PhenotypeType type = phenotype.type();
            if (type != null) {
                cancerTypeDOID = type.id();
            }
        }
        String level = entry.association().evidenceLabel();
        String direction = entry.association().responseType();

        if (drugs != null && cancerTypeString != null && cancerTypeDOID != null && level != null && direction != null) {
            return ImmutableActionableEvidence.builder()
                    .drugs(drugs)
                    .cancerTypeString(cancerTypeString)
                    .cancerTypeDOID(cancerTypeDOID)
                    .level(level)
                    .direction(direction)
                    .build();
        } else {
            return null;
        }
    }
}
