package com.hartwig.hmftools.knowledgebasegenerator.hotspot;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.knowledgebasegenerator.RefGenomeVersion;
import com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase.Source;
import com.hartwig.hmftools.knowledgebasegenerator.transvar.Transvar;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class HotspotExtractor {

    private static final Logger LOGGER = LogManager.getLogger(HotspotExtractor.class);

    private static final Set<String> ONCOKB_VALID_BIOMARKER_TYPES =
            Sets.newHashSet("missense_variant", "inframe_deletion", "inframe_insertion");

    @NotNull
    private final Transvar transvar;

    public static HotspotExtractor withRefGenome(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile)
            throws FileNotFoundException {
        return new HotspotExtractor(Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile));
    }

    private HotspotExtractor(@NotNull Transvar transvar) {
        this.transvar = transvar;
    }

    @NotNull
    public List<VariantHotspot> extractHotspots(@NotNull ViccEntry viccEntry) throws IOException, InterruptedException {
        List<VariantHotspot> allHotspots = Lists.newArrayList();
        if (Source.sourceFromKnowledgebase(viccEntry.source()) == Source.ONCOKB) {
            for (Feature feature : viccEntry.features()) {
                if (ONCOKB_VALID_BIOMARKER_TYPES.contains(feature.biomarkerType()) || isProteinAnnotation(feature.name())) {
                    List<VariantHotspot> hotspots = transvar.extractHotspotsFromProteinAnnotation(feature.geneSymbol(), feature.name());
                    LOGGER.info("Converted '{}' on gene '{}' to {} hotspot(s)", feature.name(), feature.geneSymbol(), hotspots.size());
                    allHotspots.addAll(hotspots);
                }
            }
        }

        return allHotspots;
    }

    @VisibleForTesting
    static boolean isProteinAnnotation(@NotNull String featureName) {
        String featureToTest;
        // Features could be ranges such as E102_I103del
        if (featureName.contains("_")) {
            featureToTest = featureName.split("_")[1];
            // We only want to explicitly filter on ranges that involve insertions, deletions or duplication.
            if (!(featureToTest.contains("ins") || featureToTest.contains("del") || featureToTest.contains("dup"))) {
                return false;
            }
        } else {
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
