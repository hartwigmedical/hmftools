package com.hartwig.hmftools.serve.extraction.hotspot.tools;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;
import com.hartwig.hmftools.serve.extraction.util.VCFWriterFactory;
import com.hartwig.hmftools.serve.util.AminoAcidFunctions;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class AnnotatedHotspotVCFChecker {

    private static final Logger LOGGER = LogManager.getLogger(AnnotatedHotspotVCFChecker.class);
    private static final boolean LOG_DEBUG = false;

    private static final String NO_INPUT_PROTEIN = "-";

    private final Set<String> curatedTranscripts = Sets.newHashSet();
    private final Map<String, Set<String>> annotationsRequestedForMappingPerTranscript = Maps.newHashMap();
    private final Map<String, Set<String>> annotationsRequestedForMappingPerGene = Maps.newHashMap();

    public static void main(String[] args) throws IOException {
        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String annotatedHotspotVcf = System.getProperty("user.home") + "/hmf/tmp/annotatedHotspots.vcf";
        new AnnotatedHotspotVCFChecker().run(annotatedHotspotVcf);
    }

    public void run(@NotNull String annotatedVcfFilePath) throws IOException {
        int totalCount = 0;
        int matchCount = 0;
        int whitelistedMatchCount = 0;
        int diffCount = 0;

        LOGGER.info("Loading hotspots from '{}'", annotatedVcfFilePath);
        AbstractFeatureReader<VariantContext, LineIterator> reader =
                AbstractFeatureReader.getFeatureReader(annotatedVcfFilePath, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            totalCount++;

            String[] inputParts = variant.getAttributeAsString(VCFWriterFactory.INPUT_FIELD, Strings.EMPTY).split("\\|");
            String inputGene = inputParts[0];
            String inputTranscript = inputParts[1].equals("null") ? null : inputParts[1];
            String inputProteinAnnotation = inputParts[2];

            String formattedHotspot = formatHotspot(variant);
            if (inputProteinAnnotation.equals(NO_INPUT_PROTEIN)) {
                LOGGER.debug("Skipping non-coding hotspot on '{}'", formattedHotspot);
                matchCount++;
            } else {
                List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(variant);
                MatchType match = determineMatch(inputGene, inputTranscript, inputProteinAnnotation, annotations, formattedHotspot);
                switch (match) {
                    case WHITE_LIST: {
                        whitelistedMatchCount++;
                        matchCount++;
                        break;
                    }
                    case IDENTICAL: {
                        matchCount++;
                        break;
                    }
                    case NO_MATCH: {
                        diffCount++;
                        break;
                    }
                }
            }
        }

        LOGGER.info("Done comparing {} records: {} matches (of which {} through whitelisting) and {} differences found.",
                totalCount,
                matchCount,
                whitelistedMatchCount,
                diffCount);

        checkForUnusedMappings();
    }

    @NotNull
    private static String formatHotspot(@NotNull VariantContext variant) {
        return variant.getContig() + ":" + variant.getStart() + " " + variant.getReference().getBaseString() + ">"
                + variant.getAlternateAllele(0).getBaseString();
    }

    @NotNull
    private MatchType determineMatch(@NotNull String inputGene, @Nullable String inputTranscript, @NotNull String inputProteinAnnotation,
            @NotNull List<SnpEffAnnotation> annotations, @NotNull String formattedHotspot) {
        SnpEffAnnotation specificAnnotation = annotationForTranscript(annotations, inputTranscript);
        if (specificAnnotation != null) {
            return matchOnSpecificAnnotation(inputGene, inputTranscript, inputProteinAnnotation, specificAnnotation, formattedHotspot);
        } else {
            // In case input transcript is missing or can't be found, we try to match against any transcript.
            // This could be tricky in case a variant was generated from 37 and is now being evaluated on 38 with different transcript IDs.
            return matchOnAnyTranscript(inputGene, inputProteinAnnotation, annotations);
        }
    }

    @NotNull
    private MatchType matchOnSpecificAnnotation(@NotNull String inputGene, @NotNull String inputTranscript,
            @NotNull String inputProteinAnnotation, @NotNull SnpEffAnnotation specificAnnotation, @NotNull String formattedHotspot) {
        String snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(specificAnnotation.hgvsProtein());
        if (!isSameAnnotation(inputGene, inputTranscript, inputProteinAnnotation, snpeffProteinAnnotation)) {
            LOGGER.warn("Difference on gene '{}-{}' - {} : SERVE input protein '{}' vs SnpEff protein '{}'",
                    inputGene,
                    inputTranscript,
                    formattedHotspot,
                    inputProteinAnnotation,
                    snpeffProteinAnnotation);
            return MatchType.NO_MATCH;
        } else {
            if (snpeffProteinAnnotation.equals(inputProteinAnnotation)) {
                LOGGER.debug("Identical match found on {} for '{}'", inputGene, inputProteinAnnotation);
                return MatchType.IDENTICAL;
            } else {
                LOGGER.debug("Match found on {}. '{}' and '{}' are considered identical",
                        inputGene,
                        inputProteinAnnotation,
                        snpeffProteinAnnotation);
                return MatchType.WHITE_LIST;
            }
        }
    }

    @NotNull
    private MatchType matchOnAnyTranscript(@NotNull String inputGene, @NotNull String inputProteinAnnotation,
            @NotNull List<SnpEffAnnotation> annotations) {
        boolean matchFound = false;
        String matchedSnpeffProteinAnnotation = null;
        for (SnpEffAnnotation annotation : annotations) {
            // We only want to consider transcript features with coding impact.
            if (annotation.isTranscriptFeature() && !annotation.hgvsProtein().isEmpty()) {
                String snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());
                if (isSameAnnotation(inputGene, annotation.transcript(), inputProteinAnnotation, snpeffProteinAnnotation)) {
                    matchFound = true;
                    matchedSnpeffProteinAnnotation = snpeffProteinAnnotation;
                }
            }
        }

        if (matchFound) {
            if (inputProteinAnnotation.equals(matchedSnpeffProteinAnnotation)) {
                LOGGER.debug("Found a match amongst candidate transcripts for '{}' on '{}", inputProteinAnnotation, inputGene);
                return MatchType.IDENTICAL;
            } else {
                LOGGER.debug("Match found on {}. '{}' and '{}' are considered identical",
                        inputGene,
                        inputProteinAnnotation,
                        matchedSnpeffProteinAnnotation);
                return MatchType.WHITE_LIST;
            }
        } else {
            LOGGER.warn("Could not find a match amongst candidate transcripts for '{}' on '{}'", inputProteinAnnotation, inputGene);
            return MatchType.NO_MATCH;
        }
    }

    @Nullable
    private static SnpEffAnnotation annotationForTranscript(@NotNull List<SnpEffAnnotation> annotations, @Nullable String transcript) {
        for (SnpEffAnnotation annotation : annotations) {
            if (annotation.isTranscriptFeature() && annotation.transcript().equals(transcript)) {
                return annotation;
            }
        }
        return null;
    }

    private boolean isSameAnnotation(@NotNull String gene, @NotNull String transcript, @NotNull String inputAnnotation,
            @NotNull String snpeffAnnotation) {
        String curatedInputAnnotation = curateStartCodonAnnotation(inputAnnotation);
        return curatedInputAnnotation.equals(snpeffAnnotation) || transcriptWhitelist(transcript, inputAnnotation, snpeffAnnotation)
                || geneWhitelist(gene, inputAnnotation, snpeffAnnotation) || retiredTranscriptCheck(transcript, snpeffAnnotation)
                || changedTranscriptCheck(transcript, inputAnnotation, snpeffAnnotation);
    }

    private boolean transcriptWhitelist(@NotNull String transcript, @NotNull String inputAnnotation, @NotNull String snpeffAnnotation) {
        Map<String, List<String>> transcriptMapping =
                AnnotatedHotspotCurationFactory.SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT.get(transcript);
        if (transcriptMapping != null) {
            Set<String> requestedAnnotations = annotationsRequestedForMappingPerTranscript.get(transcript);
            if (requestedAnnotations == null) {
                requestedAnnotations = Sets.newHashSet(inputAnnotation);
            } else {
                requestedAnnotations.add(inputAnnotation);
            }
            annotationsRequestedForMappingPerTranscript.put(transcript, requestedAnnotations);
            List<String> mappedAnnotations = transcriptMapping.get(inputAnnotation);
            if (mappedAnnotations != null) {
                return mappedAnnotations.contains(snpeffAnnotation);
            }
        }

        return false;
    }

    private boolean geneWhitelist(@NotNull String gene, @NotNull String inputAnnotation, @NotNull String snpeffAnnotation) {
        Map<String, List<String>> geneMapping = AnnotatedHotspotCurationFactory.SERVE_TO_SNPEFF_MAPPINGS_PER_GENE.get(gene);
        if (geneMapping != null) {
            Set<String> requestedAnnotations = annotationsRequestedForMappingPerGene.get(gene);
            if (requestedAnnotations == null) {
                requestedAnnotations = Sets.newHashSet(inputAnnotation);
            } else {
                requestedAnnotations.add(inputAnnotation);
            }
            annotationsRequestedForMappingPerGene.put(gene, requestedAnnotations);
            List<String> mappedAnnotations = geneMapping.get(inputAnnotation);
            if (mappedAnnotations != null) {
                return mappedAnnotations.contains(snpeffAnnotation);
            }
        }

        return false;
    }

    private boolean retiredTranscriptCheck(@NotNull String transcript, @NotNull String snpeffAnnotation) {
        if (AnnotatedHotspotCurationFactory.RETIRED_TRANSCRIPTS.contains(transcript) && snpeffAnnotation.isEmpty()) {
            // In case we know a transcript has been retired from coding duty in certain ref genomes we accept the diff when empty.
            curatedTranscripts.add(transcript);
            return true;
        }

        return false;
    }

    private boolean changedTranscriptCheck(@NotNull String transcript, @NotNull String inputAnnotation, @NotNull String snpeffAnnotation) {
        if (AnnotatedHotspotCurationFactory.CHANGED_TRANSCRIPTS.contains(transcript)) {
            // In case transcripts have different versions across ref genomes we assume the AA change is the same, just a different position
            // Eg p.PxxS should match for any xx
            curatedTranscripts.add(transcript);
            return inputAnnotation.substring(0, 3).equals(snpeffAnnotation.substring(0, 3))
                    && inputAnnotation.substring(inputAnnotation.length()).equals(snpeffAnnotation.substring(snpeffAnnotation.length()));
        }

        return false;
    }

    private void checkForUnusedMappings() {
        int unusedMappingCount = 0;
        Set<Map.Entry<String, Map<String, List<String>>>> transcriptEntries =
                AnnotatedHotspotCurationFactory.SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT.entrySet();
        for (Map.Entry<String, Map<String, List<String>>> entry : transcriptEntries) {
            Set<String> requestedAnnotations = annotationsRequestedForMappingPerTranscript.get(entry.getKey());
            if (requestedAnnotations == null) {
                LOGGER.warn("No annotation mapping requested at all for '{}'", entry.getKey());
                unusedMappingCount += entry.getValue().keySet().size();
            } else {
                for (String annotationKey : entry.getValue().keySet()) {
                    if (!requestedAnnotations.contains(annotationKey)) {
                        LOGGER.warn("Unused annotation configured for '{}': '{}'", entry.getKey(), annotationKey);
                        unusedMappingCount++;
                    }
                }
            }
        }

        LOGGER.info("Analyzed usage of mapping configuration. Found {} unused mappings.", unusedMappingCount);

        int unusedTranscriptCurationCount = 0;
        for (String transcript : AnnotatedHotspotCurationFactory.RETIRED_TRANSCRIPTS) {
            if (!curatedTranscripts.contains(transcript)) {
                LOGGER.warn("Transcript '{}' labeled as 'retired' has not been used in curation", transcript);
                unusedTranscriptCurationCount++;
            }
        }

        for (String transcript : AnnotatedHotspotCurationFactory.CHANGED_TRANSCRIPTS) {
            if (!curatedTranscripts.contains(transcript)) {
                LOGGER.warn("Transcript '{}' labeled as 'changed' has not been used in curation", transcript);
                unusedTranscriptCurationCount++;
            }
        }

        LOGGER.info("Analyzed usage of transcript curation. Found {} unused curation keys.", unusedTranscriptCurationCount);
    }

    @NotNull
    private static String curateStartCodonAnnotation(@NotNull String serveAnnotation) {
        if (serveAnnotation.startsWith("p.M1") && serveAnnotation.length() == 5) {
            return "p.M1?";
        } else {
            return serveAnnotation;
        }
    }

    private enum MatchType {
        IDENTICAL,
        WHITE_LIST,
        NO_MATCH
    }
}
