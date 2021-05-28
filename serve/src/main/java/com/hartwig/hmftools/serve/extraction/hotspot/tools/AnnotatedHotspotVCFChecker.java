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

    private static final Set<String> RETIRED_TRANSCRIPTS = createRetiredTranscriptWhitelist();
    private static final Set<String> CHANGED_TRANSCRIPTS = createChangedTranscriptWhitelist();
    private static final Map<String, Map<String, List<String>>> SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT = createMappings();

    private final Set<String> curatedTranscripts = Sets.newHashSet();
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
            String snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(specificAnnotation.hgvsProtein());
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
            // In case input transcript is missing or can't be found, we try to match against any transcript.
            // This could be tricky in case a variant was generated from 37 and is now being evaluated on 38 with different transcript IDs.
            boolean matchFound = false;
            for (SnpEffAnnotation annotation : annotations) {
                // We only want to consider transcript features with coding impact.
                if (annotation.isTranscriptFeature() && !annotation.hgvsProtein().isEmpty()) {
                    String snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());
                    if (isSameAnnotation(annotation.transcript(), inputProteinAnnotation, snpeffProteinAnnotation)) {
                        matchFound = true;
                    }
                }
            }

            if (matchFound) {
                LOGGER.debug("Found a match amongst candidate transcripts for '{}' on '{}", inputProteinAnnotation, inputGene);
                return MatchType.IDENTICAL;
            } else if (inputTranscript != null && RETIRED_TRANSCRIPTS.contains(inputTranscript)) {
                LOGGER.debug("Transcript '{}' has retired. Whitelisting entry '{}' on '{}'",
                        inputTranscript,
                        inputProteinAnnotation,
                        inputGene);
                curatedTranscripts.add(inputTranscript);
                return MatchType.WHITE_LIST;
            } else {
                LOGGER.warn("Could not find a match amongst candidate transcripts for '{}' on '{}'", inputProteinAnnotation, inputGene);
                return MatchType.NO_MATCH;
            }
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

        if (RETIRED_TRANSCRIPTS.contains(transcript) && snpeffAnnotation.isEmpty()) {
            // In case we know a transcript has been retired from coding duty in certain ref genomes we accept the diff when empty.
            curatedTranscripts.add(transcript);
            return true;
        }

        if (CHANGED_TRANSCRIPTS.contains(transcript)) {
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

        int unusedTranscriptCurationCount = 0;
        for (String transcript : RETIRED_TRANSCRIPTS) {
            if (!curatedTranscripts.contains(transcript)) {
                LOGGER.warn("Transcript '{}' labeled as 'retired' has not been used in curation", transcript);
                unusedTranscriptCurationCount++;
            }
        }

        for (String transcript : CHANGED_TRANSCRIPTS) {
            if (!curatedTranscripts.contains(transcript)) {
                LOGGER.warn("Transcript '{}' labeled as 'changed' has not been used in curation", transcript);
                unusedTranscriptCurationCount++;
            }
        }

        LOGGER.info("Analyzed usage of transcript curation. Found {} unused curation keys.", unusedTranscriptCurationCount);
    }

    @NotNull
    private static Set<String> createRetiredTranscriptWhitelist() {
        Set<String> transcripts = Sets.newHashSet();

        // This transcript on RHOA is coding in v37 but non-coding in v38
        transcripts.add("ENST00000431929");

        // This transcript on CDKN2A has retired in v38
        transcripts.add("ENST00000361570");

        // This transcript on BCL2L12 has completely changed in 38 and can't match anymore.
        transcripts.add("ENST00000246784");

        return transcripts;
    }

    @NotNull
    private static Set<String> createChangedTranscriptWhitelist() {
        Set<String> transcripts = Sets.newHashSet();

        // This MYC transcript is v2 in 37 and v6 in 38 and they differ by 15 AA.
        transcripts.add("ENST00000377970");

        // This TCF7L2 transcript in v1 in 37 and v5 in 38 and they differ by 5 AA.
        transcripts.add("ENST00000543371");

        // This FGFR2 transcript is v6 in 37 and v10 in 38 and they differ by 36 AA (spread across transcript).
        transcripts.add("ENST00000351936");

        // This MED12 transcript is v6 in 37 and v10 in 38 and they differ with 53 AA (spread across transcript)
        transcripts.add("ENST00000333646");

        return transcripts;
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
