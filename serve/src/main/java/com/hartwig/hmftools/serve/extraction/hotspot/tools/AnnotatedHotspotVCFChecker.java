package com.hartwig.hmftools.serve.extraction.hotspot.tools;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;
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

    private static final Map<String, Map<String, List<String>>> SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT = createMappings();

    private final Map<String, Set<String>> annotationsRequestedForMappingPerTranscript = Maps.newHashMap();

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

            String[] inputParts = variant.getAttributeAsString("input", Strings.EMPTY).split("\\|");
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
        if (inputTranscript != null) {
            SnpEffAnnotation annotation = annotationForTranscript(annotations, inputTranscript);

            if (annotation != null) {
                String snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());
                if (!isSameAnnotation(inputTranscript, inputProteinAnnotation, snpeffProteinAnnotation)) {
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
            } else {
                LOGGER.warn("Could not find snpeff annotation for '{}' on '{}'!", inputTranscript, inputGene);
                return MatchType.NO_MATCH;
            }
        } else {
            // In case input transcript is missing we try to match against any transcript.
            boolean matchFound = false;
            for (SnpEffAnnotation annotation : annotations) {
                if (annotation.isTranscriptFeature()) {
                    String snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());
                    if (isSameAnnotation(annotation.transcript(), inputProteinAnnotation, snpeffProteinAnnotation)) {
                        matchFound = true;
                    }
                }
            }

            if (matchFound) {
                LOGGER.debug("Found a match amongst candidate transcripts for '{}' on '{}", inputProteinAnnotation, inputGene);
                return MatchType.IDENTICAL;
            } else {
                LOGGER.warn("Could not find a match amongst candidate transcripts for '{}' on '{}'", inputProteinAnnotation, inputGene);
                return MatchType.NO_MATCH;
            }
        }
    }

    @Nullable
    private static SnpEffAnnotation annotationForTranscript(@NotNull List<SnpEffAnnotation> annotations, @NotNull String transcript) {
        for (SnpEffAnnotation annotation : annotations) {
            if (annotation.isTranscriptFeature() && annotation.transcript().equals(transcript)) {
                return annotation;
            }
        }
        return null;
    }

    private boolean isSameAnnotation(@NotNull String transcript, @NotNull String inputAnnotation, @NotNull String snpeffAnnotation) {
        String curatedInputAnnotation = curateStartCodonAnnotation(inputAnnotation);
        if (curatedInputAnnotation.equals(snpeffAnnotation)) {
            return true;
        }

        Map<String, List<String>> transcriptMapping = SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT.get(transcript);
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

    private void checkForUnusedMappings() {
        int unusedMappingCount = 0;
        for (Map.Entry<String, Map<String, List<String>>> entry : SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT.entrySet()) {
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
    }

    @NotNull
    private static Map<String, Map<String, List<String>>> createMappings() {
        Map<String, Map<String, List<String>>> serveToSnpEffMappings = Maps.newHashMap();

        serveToSnpEffMappings.put("ENST00000288135", createKITMap());
        serveToSnpEffMappings.put("ENST00000275493", createEGFRMap());
        serveToSnpEffMappings.put("ENST00000269571", createERBB2Map());
        serveToSnpEffMappings.put("ENST00000263967", createPIK3CAMap());
        serveToSnpEffMappings.put("ENST00000374690", createARMap());
        serveToSnpEffMappings.put("ENST00000231790", createMLH1Map());

        return serveToSnpEffMappings;
    }

    @NotNull
    private static Map<String, List<String>> createKITMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.S501_A502insAY", Lists.newArrayList("p.A502_Y503insYA"));
        map.put("p.V560del", Lists.newArrayList("p.V559del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createEGFRMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.V769_D770insASV", Lists.newArrayList("p.A767_V769dup"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createERBB2Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.A771_Y772insYVMA", Lists.newArrayList("p.Y772_V773insVMAY"));
        map.put("p.Y772_A775dup", Lists.newArrayList("p.Y772_V773insVMAY"));
        map.put("p.M774_A775insAYVM", Lists.newArrayList("p.A775_G776insYVMA"));
        map.put("p.G778_P780dup", Lists.newArrayList("p.G778_S779insSPG"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createPIK3CAMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.E109del", Lists.newArrayList("p.E110del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createARMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.Q86del", Lists.newArrayList("p.Q87del", "p.Q88del", "p.Q89del", "p.Q90del", "p.Q91del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createMLH1Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.K618del", Lists.newArrayList("p.K616del", "p.K617del"));
        return map;
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
