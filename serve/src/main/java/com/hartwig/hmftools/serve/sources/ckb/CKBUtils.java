package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.extraction.signature.SignatureName;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CKBUtils {

    private static final Logger LOGGER = LogManager.getLogger(CKBUtils.class);

    public static void printExtractionResults(@NotNull Map<CkbEntry, ExtractResult> resultsPerEntry) {
        List<Variant> featuresWithoutGenomicEvents = Lists.newArrayList();

        int totalFeatureCount = 0;
        int featuresWithHotspotsCount = 0;
        int totalHotspotsCount = 0;
        int featuresWithCodonCount = 0;
        int totalCodonCount = 0;
        int featuresWithExonCount = 0;
        int totalExonCount = 0;

        int featuresWithGeneLevelEventCount = 0;
        int featuresWithCopyNumberCount = 0;
        int featuresWithFusionCount = 0;
        int featuresWithSignatureCount = 0;

        for (Map.Entry<CkbEntry, ExtractResult> resultPerEntry : resultsPerEntry.entrySet()) {
            CkbEntry entry = resultPerEntry.getKey();
            ExtractResult result = resultPerEntry.getValue();
            for (Variant variant : entry.variants()) {
                List<VariantHotspot> hotspotsForFeature = result.hotspotsPerFeature().get(variant);
                List<CodonAnnotation> codonsForFeature = result.codonsPerFeature().get(variant);
                List<ExonAnnotation> exonsForFeature = result.exonsPerFeature().get(variant);
                GeneLevelAnnotation geneLevelEventForFeature = result.geneLevelEventsPerFeature().get(variant);
                KnownCopyNumber ampDelForFeature = result.ampsDelsPerFeature().get(variant);
                KnownFusionPair fusionForFeature = result.fusionsPerFeature().get(variant);
                SignatureName signatureForFeature = result.signaturesPerFeature().get(variant);

                if (hotspotsForFeature == null && codonsForFeature == null && exonsForFeature == null && geneLevelEventForFeature == null
                        && ampDelForFeature == null && fusionForFeature == null && signatureForFeature == null) {
                    //                            if (feature.type() != EventType.COMBINED && feature.type() != EventType.COMPLEX) {
                    //                                // For both combined and complex events we expect no genomic events to be derived.
                    //                                featuresWithoutGenomicEvents.add(variant);
                    //                            }
                } else {
                    if (hotspotsForFeature != null) {
                        featuresWithHotspotsCount++;
                        totalHotspotsCount += hotspotsForFeature.size();
                    }

                    if (codonsForFeature != null) {
                        featuresWithCodonCount++;
                        totalCodonCount += codonsForFeature.size();
                    }

                    if (exonsForFeature != null) {
                        featuresWithExonCount++;
                        totalExonCount += exonsForFeature.size();
                    }

                    if (geneLevelEventForFeature != null) {
                        featuresWithGeneLevelEventCount++;
                    }

                    if (ampDelForFeature != null) {
                        featuresWithCopyNumberCount++;
                    }

                    if (fusionForFeature != null) {
                        featuresWithFusionCount++;
                    }

                    if (signatureForFeature != null) {
                        featuresWithSignatureCount++;
                    }
                }

                totalFeatureCount++;
            }
        }

        if (!featuresWithoutGenomicEvents.isEmpty()) {
            LOGGER.warn("No genomic events derived for {} features!", featuresWithoutGenomicEvents.size());
            for (Variant variant : featuresWithoutGenomicEvents) {
                LOGGER.debug(" No genomic events derived from '{}' in '{}'", variant.variant(), variant.gene().geneSymbol());
            }
        }

        LOGGER.info("Analysis performed on {} features in {} CKB entries", totalFeatureCount, resultsPerEntry.size());
        LOGGER.info(" Extracted {} hotspots for {} features", totalHotspotsCount, featuresWithHotspotsCount);
        LOGGER.info(" Extracted {} codons for {} features", totalCodonCount, featuresWithCodonCount);
        LOGGER.info(" Extracted {} exons for {} features", totalExonCount, featuresWithExonCount);
        LOGGER.info(" Extracted {} gene level events", featuresWithGeneLevelEventCount);
        LOGGER.info(" Extracted {} known amps and dels", featuresWithCopyNumberCount);
        LOGGER.info(" Extracted {} known fusions pairs", featuresWithFusionCount);
        LOGGER.info(" Extracted {} signatures", featuresWithSignatureCount);
    }
}

