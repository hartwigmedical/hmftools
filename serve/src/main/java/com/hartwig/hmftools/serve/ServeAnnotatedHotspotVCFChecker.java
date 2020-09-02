package com.hartwig.hmftools.serve;

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

public class ServeAnnotatedHotspotVCFChecker {

    private static final Logger LOGGER = LogManager.getLogger(ServeAnnotatedHotspotVCFChecker.class);
    private static final boolean LOG_DEBUG = false;

    private static final String NO_INPUT_PROTEIN = "-";

    private static final Map<String, Map<String, List<String>>> SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT = createMappings();

    private final Map<String, Set<String>> annotationsRequestedForMappingPerTranscript = Maps.newHashMap();

    public static void main(String[] args) throws IOException {
        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String annotatedHotspotVcf = System.getProperty("user.home") + "/hmf/tmp/annotatedHotspots.vcf";
        new ServeAnnotatedHotspotVCFChecker().run(annotatedHotspotVcf);
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

        serveToSnpEffMappings.put("ENST00000257290", createPDGFRAMap());
        serveToSnpEffMappings.put("ENST00000288135", createKITMap());
        serveToSnpEffMappings.put("ENST00000274335", createPIK3R1Map());
        serveToSnpEffMappings.put("ENST00000275493", createEGFRMap());
        serveToSnpEffMappings.put("ENST00000288602", createBRAFMap());
        serveToSnpEffMappings.put("ENST00000277541", createNOTCH1Map());
        serveToSnpEffMappings.put("ENST00000358487", createFGFR2Map());
        serveToSnpEffMappings.put("ENST00000241453", createFLT3Map());
        serveToSnpEffMappings.put("ENST00000262367", createCREBBPMap());
        serveToSnpEffMappings.put("ENST00000269571", createERBB2Map());
        serveToSnpEffMappings.put("ENST00000357654", createBRCA1Map());
        serveToSnpEffMappings.put("ENST00000377045", createARAFMap());
        serveToSnpEffMappings.put("ENST00000460911", createEZH2Map()); // Note: This is not the canonical transcript for EZH2
        serveToSnpEffMappings.put("ENST00000263967", createPIK3CAMap());
        serveToSnpEffMappings.put("ENST00000374690", createARMap());

        return serveToSnpEffMappings;
    }

    @NotNull
    private static Map<String, List<String>> createPDGFRAMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.R841_D842delinsKI", Lists.newArrayList("p.RD841KI"));
        map.put("p.D842_I843delinsIM", Lists.newArrayList("p.DI842IM"));
        map.put("p.D842_H845del", Lists.newArrayList("p.I843_D846del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createKITMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.S501_A502insAY", Lists.newArrayList("p.A502_Y503insYA"));
        map.put("p.K550_W557del", Lists.newArrayList("p.P551_K558del"));
        map.put("p.V555_V559del", Lists.newArrayList("p.Q556_V560del"));
        map.put("p.P577_W582delinsPYD", Lists.newArrayList("p.H580_W582del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createPIK3R1Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.I559_D560insDKRMNS", Lists.newArrayList("p.K561_R562insRMNSDK"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createEGFRMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.K745_E746insIPVAIK", Lists.newArrayList("p.I740_K745dup"));
        map.put("p.S768_V769insVAS", Lists.newArrayList("p.V769_D770insASV"));
        map.put("p.D770_N771insD", Lists.newArrayList("p.D770dup"));
        map.put("p.V769_D770insASV", Lists.newArrayList("p.A767_V769dup"));
        map.put("p.D770_N771insSVD", Lists.newArrayList("p.S768_D770dup"));
        map.put("p.D770_N771insNPG", Lists.newArrayList("p.P772_H773insGNP"));
        map.put("p.H773_V774insH", Lists.newArrayList("p.H773dup"));
        map.put("p.H773_V774insNPH", Lists.newArrayList("p.N771_H773dup"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createBRAFMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.T599_V600insV", Lists.newArrayList("p.V600dup"));
        map.put("p.T599_V600insEAT", Lists.newArrayList("p.A598_T599insTEA"));
        map.put("p.T599_V600insETT", Lists.newArrayList("p.A598_T599insTET"));
        map.put("p.L485_Q494del", Lists.newArrayList("p.N486_L495del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createNOTCH1Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.P2415del", Lists.newArrayList("p.P2414del", "p.P2413del", "p.P2412del", "p.P2411del"));
        map.put("p.V1578del", Lists.newArrayList("p.V1577del", "p.V1576del", "p.V1575del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createFGFR2Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.S267_D273dup", Lists.newArrayList("p.D273_V274insSTVVGGD"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createFLT3Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.W603_E604insDREYEYDLKW", Lists.newArrayList("p.K602_W603insWDREYEYDLK"));
        map.put("p.Y599_D600insEYEYEYEY", Lists.newArrayList("p.E598_Y599insYEYEYEYE"));
        map.put("p.Y599_D600insGLYVDFREYEY", Lists.newArrayList("p.E598_Y599insYGLYVDFREYE"));
        map.put("p.Y599_D600insSTDNEYFYVDFREYEY", Lists.newArrayList("p.E598_Y599insYSTDNEYFYVDFREYE"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createCREBBPMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.S1680del", Lists.newArrayList("p.S1679del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createERBB2Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.A771_Y772insYVMA", Lists.newArrayList("p.Y772_V773insVMAY"));
        map.put("p.M774_A775insAYVM", Lists.newArrayList("p.A775_G776insYVMA"));
        map.put("p.A775_G776insYVMA", Lists.newArrayList("p.Y772_A775dup"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createBRCA1Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.V1688del", Lists.newArrayList("p.V1687del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createARAFMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.Q347_A348del", Lists.newArrayList("p.Q349_A350del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createEZH2Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.T678_R679delinsKK", Lists.newArrayList("p.TR678KK"));
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
