package com.hartwig.hmftools.serve.vicc.hotspot;

import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_DELETION;
import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_DUPLICATION;
import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_FRAMESHIFT_SUFFIX;
import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_FRAMESHIFT_SUFFIX_WITH_STOP_GAINED;
import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_INSERTION;
import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_RANGE_INDICATOR;

import java.io.FileNotFoundException;
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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HotspotExtractor {

    private static final Logger LOGGER = LogManager.getLogger(HotspotExtractor.class);

    private static final int MAX_INFRAME_BASE_LENGTH = 50;

    @NotNull
    private final ProteinResolver proteinResolver;

    @NotNull
    public static HotspotExtractor transvarWithRefGenome(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile)
            throws FileNotFoundException {
        LOGGER.info("Creating hotspot extractor with ref genome version '{}' and fasta path '{}'", refGenomeVersion, refGenomeFastaFile);
        return new HotspotExtractor(Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile));
    }

    @NotNull
    public static HotspotExtractor dummy() {
        return new HotspotExtractor(new ProteinResolver() {
            @NotNull
            @Override
            public List<VariantHotspot> extractHotspotsFromProteinAnnotation(@NotNull final String gene,
                    @Nullable final String specificTranscript, @NotNull final String proteinAnnotation) {
                return Lists.newArrayList();
            }

            @NotNull
            @Override
            public Set<String> unresolvedProteinAnnotations() {
                return Sets.newHashSet();
            }
        });
    }

    private HotspotExtractor(@NotNull final ProteinResolver proteinResolver) {
        this.proteinResolver = proteinResolver;
    }

    @NotNull
    public Map<Feature, List<VariantHotspot>> extractHotspots(@NotNull ViccEntry viccEntry) {
        Map<Feature, List<VariantHotspot>> allHotspotsPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            if (isResolvableProteinAnnotation(feature.proteinAnnotation())) {
                List<VariantHotspot> hotspots = proteinResolver.extractHotspotsFromProteinAnnotation(feature.geneSymbol(),
                        viccEntry.transcriptId(),
                        feature.proteinAnnotation());

                allHotspotsPerFeature.put(feature, hotspots);
            }
        }

        return allHotspotsPerFeature;
    }

    @NotNull
    public Set<String> unresolvedProteinAnnotations() {
        return proteinResolver.unresolvedProteinAnnotations();
    }

    @VisibleForTesting
    static boolean isResolvableProteinAnnotation(@NotNull String feature) {
        try {
            if (isFrameshift(feature)) {
                return false;
            } else if (feature.contains(HGVS_RANGE_INDICATOR)) {
                return isValidRangeMutation(feature);
            } else {
                return isValidSingleCodonMutation(feature);
            }
        } catch (Exception exception) {
            LOGGER.warn("Could not determine whether feature is protein annotation due to '{}'", exception.getMessage(), exception);
            return false;
        }
    }

    private static boolean isFrameshift(@NotNull String feature) {
        return feature.endsWith(HGVS_FRAMESHIFT_SUFFIX) || feature.endsWith(HGVS_FRAMESHIFT_SUFFIX_WITH_STOP_GAINED);
    }

    private static boolean isValidRangeMutation(@NotNull String feature) {
        assert feature.contains(HGVS_RANGE_INDICATOR);

        // Features could be ranges such as E102_I103del. We whitelist specific feature types when analyzing a range.
        String[] featureParts = feature.split(HGVS_RANGE_INDICATOR);
        String featureStartPart = featureParts[0];
        String featureEndPart = featureParts[1];
        if (featureEndPart.contains(HGVS_INSERTION) || featureEndPart.contains(HGVS_DUPLICATION)
                || featureEndPart.contains(HGVS_DELETION)) {
            int indexOfEvent;
            if (featureEndPart.contains(HGVS_INSERTION)) {
                indexOfEvent = featureEndPart.indexOf(HGVS_INSERTION);
            } else if (featureEndPart.contains(HGVS_DUPLICATION)) {
                indexOfEvent = featureEndPart.indexOf(HGVS_DUPLICATION);
            } else {
                indexOfEvent = featureEndPart.indexOf(HGVS_DELETION);
            }

            long start = Long.parseLong(featureStartPart.substring(1));
            long end = Long.parseLong(featureEndPart.substring(1, indexOfEvent));
            return 3 * (1 + end - start) <= MAX_INFRAME_BASE_LENGTH;
        } else {
            return false;
        }
    }

    private static boolean isValidSingleCodonMutation(@NotNull String feature) {
        if (feature.contains(HGVS_INSERTION)) {
            // Insertions are only allowed in a range, since we need to know where to insert the sequence exactly.
            return false;
        }

        // Features are expected to look something like V600E (1 char - N digits - M chars)
        if (feature.length() < 3) {
            return false;
        }

        if (!Character.isLetter(feature.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(feature.charAt(1))) {
            return false;
        }

        boolean haveObservedNonDigit = !Character.isDigit(feature.charAt(2));
        int firstNotDigit = -1;
        for (int i = 3; i < feature.length(); i++) {
            char charToEvaluate = feature.charAt(i);
            if (haveObservedNonDigit && Character.isDigit(charToEvaluate)) {
                return false;
            }
            boolean isDigit = Character.isDigit(charToEvaluate);
            if (!isDigit && firstNotDigit == -1) {
                firstNotDigit = i;
            }

            haveObservedNonDigit = haveObservedNonDigit || !isDigit;
        }

        if (!haveObservedNonDigit) {
            return false;
        }

        String newAminoAcid = feature.substring(firstNotDigit);
        return !newAminoAcid.equals("X");
    }
}
