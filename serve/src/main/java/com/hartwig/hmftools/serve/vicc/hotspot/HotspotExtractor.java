package com.hartwig.hmftools.serve.vicc.hotspot;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.transvar.Transvar;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class HotspotExtractor {

    private static final Logger LOGGER = LogManager.getLogger(HotspotExtractor.class);

    private static final String FEATURE_RANGE_INDICATOR = "_";
    private static final Set<String> VALID_FEATURE_RANGES = Sets.newHashSet("ins", "dup", "del");
    private static final String FRAMESHIFT_FEATURE_SUFFIX = "fs";

    @NotNull
    private final Transvar transvar;
    private final boolean transvarEnabled; // TODO Remove -> Just for testing during development
    @NotNull
    private final Set<String> unresolvableFeatures = Sets.newHashSet();

    @NotNull
    public static HotspotExtractor withRefGenome(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile,
            boolean transvarEnabled) throws FileNotFoundException {
        LOGGER.info("Creating hotspot extractor with ref genome version '{}' and fasta path '{}'", refGenomeVersion, refGenomeFastaFile);
        return new HotspotExtractor(Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile), transvarEnabled);
    }

    private HotspotExtractor(@NotNull final Transvar transvar, final boolean transvarEnabled) {
        this.transvar = transvar;
        this.transvarEnabled = transvarEnabled;
    }

    @NotNull
    public Map<Feature, List<VariantHotspot>> extractHotspots(@NotNull ViccEntry viccEntry) throws IOException, InterruptedException {
        Map<Feature, List<VariantHotspot>> allHotspotsPerFeature = Maps.newHashMap();
        if (viccEntry.source() == ViccSource.ONCOKB) {
            for (Feature feature : viccEntry.features()) {
                String featureKey = feature.geneSymbol() + ":p." + feature.name() + " - " + viccEntry.transcriptId();
                if (isProteinAnnotation(feature.name())) {
                    List<VariantHotspot> hotspots = Lists.newArrayList();
                    if (transvarEnabled) {
                        hotspots = transvar.extractHotspotsFromProteinAnnotation(feature.geneSymbol(),
                                viccEntry.transcriptId(),
                                feature.name());
                        LOGGER.debug("Converted '{}' to {} hotspot(s)", featureKey, hotspots.size());
                        if (hotspots.isEmpty()) {
                            unresolvableFeatures.add(featureKey);
                        }
                    }
                    allHotspotsPerFeature.put(feature, hotspots);
                }
            }
        }

        return allHotspotsPerFeature;
    }

    @NotNull
    public Set<String> unresolvableFeatures() {
        return unresolvableFeatures;
    }

    @VisibleForTesting
    static boolean isProteinAnnotation(@NotNull String featureName) {
        if (featureName.contains(FEATURE_RANGE_INDICATOR)) {
            // Features could be ranges such as E102_I103del. We whitelist specific feature types when analyzing a range.
            String featureToTest = featureName.split(FEATURE_RANGE_INDICATOR)[1];
            boolean validFeatureFound = false;
            for (String validFeature : VALID_FEATURE_RANGES) {
                if (featureToTest.contains(validFeature)) {
                    validFeatureFound = true;
                    break;
                }
            }
            return validFeatureFound;
        } else if (featureName.endsWith(FRAMESHIFT_FEATURE_SUFFIX)) {
            // Frameshifts are ignored for hotspot determination
            return false;
        } else {
            if (featureName.contains("ins")) {
                // "ins" is only allowed in a range, since we need to know where to insert the sequence exactly.
                return false;
            }

            // Features are expected to look something like V600E (1 char - N digits - M chars)
            if (featureName.length() < 3) {
                return false;
            }

            if (!Character.isLetter(featureName.charAt(0))) {
                return false;
            }

            if (!Character.isDigit(featureName.charAt(1))) {
                return false;
            }

            boolean haveObservedNonDigit = !Character.isDigit(featureName.charAt(2));
            for (int i = 3; i < featureName.length(); i++) {
                char charToEvaluate = featureName.charAt(i);
                if (haveObservedNonDigit && Character.isDigit(charToEvaluate)) {
                    return false;
                }
                haveObservedNonDigit = haveObservedNonDigit || !Character.isDigit(charToEvaluate);
            }

            return haveObservedNonDigit;
        }
    }
}
