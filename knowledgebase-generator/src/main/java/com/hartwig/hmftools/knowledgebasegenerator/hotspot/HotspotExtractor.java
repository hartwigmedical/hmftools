package com.hartwig.hmftools.knowledgebasegenerator.hotspot;

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

    public static HotspotExtractor withRefGenome(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile) {
        return new HotspotExtractor(new Transvar(refGenomeVersion, refGenomeFastaFile));
    }

    private HotspotExtractor(@NotNull Transvar transvar) {
        this.transvar = transvar;
    }

    @NotNull
    public List<VariantHotspot> extractHotspots(@NotNull ViccEntry viccEntry) throws IOException, InterruptedException {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        if (Source.sourceFromKnowledgebase(viccEntry.source()) == Source.ONCOKB) {
            for (Feature feature : viccEntry.features()) {
                if (ONCOKB_VALID_BIOMARKER_TYPES.contains(feature.biomarkerType())) {
//                    LOGGER.info("Converting '{}' on {} with name '{}'", feature.biomarkerType(), feature.geneSymbol(), feature.name());
                    //                    hotspots.addAll(transvar.extractHotspotsFromProteinAnnotation(feature.geneSymbol(), feature.name()));
                } else if (isProteinAnnotation(feature.name())) {
                    LOGGER.info("Attempt to convert '{}' on {}", feature.name(), feature.geneSymbol());
                } else {
                    LOGGER.info("Skipping feature interpretation of '{}' on gene '{}' with biomarker type '{}'",
                            feature.name(),
                            feature.geneSymbol(),
                            feature.biomarkerType());
                }
            }
        }

        return hotspots;
    }

    @VisibleForTesting
    static boolean isProteinAnnotation(@NotNull String featureName) {
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
