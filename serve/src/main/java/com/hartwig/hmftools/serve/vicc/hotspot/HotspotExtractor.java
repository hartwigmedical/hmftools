package com.hartwig.hmftools.serve.vicc.hotspot;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
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

    private static final Set<String> ONCOKB_VALID_BIOMARKER_TYPES =
            Sets.newHashSet("missense_variant", "inframe_deletion", "inframe_insertion");

    private static final String FEATURE_RANGE_INDICATOR = "_";
    private static final Set<String> VALID_FEATURE_RANGES = Sets.newHashSet("ins", "dup", "del");
    private static final String FRAMESHIFT_FEATURE_SUFFIX = "fs";

    @NotNull
    private final Transvar transvar;
    @NotNull
    private final Set<String> unresolvableFeatures = Sets.newHashSet();

    @NotNull
    public static HotspotExtractor withRefGenome(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile)
            throws FileNotFoundException {
        LOGGER.info("Creating hotspot extractor with ref genome version '{}' and fasta path '{}'", refGenomeVersion, refGenomeFastaFile);
        return new HotspotExtractor(Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile));
    }

    private HotspotExtractor(@NotNull Transvar transvar) {
        this.transvar = transvar;
    }

    @NotNull
    public Map<Feature, List<VariantHotspot>> extractHotspots(@NotNull ViccEntry viccEntry) throws IOException, InterruptedException {
        Map<Feature, List<VariantHotspot>> allHotspotsPerFeature = Maps.newHashMap();
        if (viccEntry.source() == ViccSource.ONCOKB) {
            for (Feature feature : viccEntry.features()) {
                String featureKey = feature.geneSymbol() + ":p." + feature.name();
                if (ONCOKB_VALID_BIOMARKER_TYPES.contains(feature.biomarkerType()) || isProteinAnnotation(feature.name())) {
                    List<VariantHotspot> hotspots = transvar.extractHotspotsFromProteinAnnotation(feature.geneSymbol(), feature.name());
                    LOGGER.info("Converted '{}' to {} hotspot(s)", featureKey, hotspots.size());
                    if (hotspots.isEmpty()) {
                        unresolvableFeatures.add(featureKey);
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
        String featureToTest;
        if (featureName.contains(FEATURE_RANGE_INDICATOR)) {
            // Features could be ranges such as E102_I103del. We whitelist specific feature types when analyzing a range.
            featureToTest = featureName.split(FEATURE_RANGE_INDICATOR)[1];
            boolean validFeatureFound = false;
            for (String validFeature : VALID_FEATURE_RANGES) {
                if (featureToTest.contains(validFeature)) {
                    validFeatureFound = true;
                    break;
                }
            }
            if (!validFeatureFound) {
                return false;
            }
        } else if (featureName.endsWith(FRAMESHIFT_FEATURE_SUFFIX)) {
            // Frameshifts are ignored for hotspot determination
            return false;
        }
        else {
            featureToTest = featureName;
        }

        // Features are expected to look something like V600E (1 char - N digits - M chars)
        if (featureToTest.length() < 3) {
            return false;
        }

        if (!Character.isLetter(featureToTest.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(featureToTest.charAt(1))) {
            return false;
        }

        boolean haveObservedNonDigit = !Character.isDigit(featureToTest.charAt(2));
        for (int i = 3; i < featureToTest.length(); i++) {
            char charToEvaluate = featureToTest.charAt(i);
            if (haveObservedNonDigit && Character.isDigit(charToEvaluate)) {
                return false;
            }
            haveObservedNonDigit = haveObservedNonDigit || !Character.isDigit(charToEvaluate);
        }

        return haveObservedNonDigit;
    }
}
