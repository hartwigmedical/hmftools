package com.hartwig.hmftools.serve;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class ServeAnnotatedHotspotVCFChecker {

    private static final Logger LOGGER = LogManager.getLogger(ServeAnnotatedHotspotVCFChecker.class);

    private static final Map<String, Map<String, List<String>>> SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT = createMappings();

    private final Map<String, Set<String>> annotationsRequestedForMappingPerTranscript = Maps.newHashMap();

    public static void main(String[] args) throws IOException {
        //        Configurator.setRootLevel(Level.DEBUG);

        String annotatedHotspotVcf = System.getProperty("user.home") + "/hmf/tmp/annotatedHotspots.vcf";
        new ServeAnnotatedHotspotVCFChecker().run(annotatedHotspotVcf);
    }

    public void run(@NotNull String annotatedVcfFilePath) throws IOException {
        int totalCount = 0;
        int matchCount = 0;
        int whitelistedMatchCount = 0;
        int diffCount = 0;

        LOGGER.info("Loading hotspots from '{}'", annotatedVcfFilePath);
        AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(annotatedVcfFilePath, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            String[] featureParts = variant.getAttributeAsString("feature", Strings.EMPTY).split("\\|");
            String featureGene = featureParts[0];
            String featureTranscript = featureParts[1].equals("null") ? null : featureParts[1];
            String featureProteinAnnotation = featureParts[2];
            if (!featureProteinAnnotation.endsWith("fs")) {
                totalCount++;

                List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(variant);

                if (featureTranscript != null) {
                    SnpEffAnnotation annotation = annotationForTranscript(annotations, featureTranscript);

                    if (annotation != null) {
                        String snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());
                        if (!isSameAnnotation(featureTranscript, featureProteinAnnotation, snpeffProteinAnnotation)) {
                            LOGGER.warn("Difference on gene '{}-{}' - {}:{} {}->{} : SERVE input protein '{}' vs SnpEff protein '{}'",
                                    featureGene,
                                    featureTranscript,
                                    variant.getContig(),
                                    variant.getStart(),
                                    variant.getReference().getBaseString(),
                                    variant.getAlternateAllele(0).getBaseString(),
                                    featureProteinAnnotation,
                                    snpeffProteinAnnotation);
                            diffCount++;
                        } else {
                            if (snpeffProteinAnnotation.equals(featureProteinAnnotation)) {
                                LOGGER.debug("Identical match found on {} for '{}'", featureGene, featureProteinAnnotation);
                            } else {
                                whitelistedMatchCount++;
                                LOGGER.debug("Match found on {}. '{}' and '{}' are considered identical",
                                        featureGene,
                                        featureProteinAnnotation,
                                        snpeffProteinAnnotation);
                            }
                            matchCount++;
                        }
                    } else {
                        LOGGER.warn("Could not find snpeff annotation for '{}' on '{}'!", featureTranscript, featureGene);
                        diffCount++;
                    }
                } else {
                    boolean matchFound = false;
                    for (SnpEffAnnotation annotation : annotations) {
                        if (annotation.isTranscriptFeature()) {
                            String snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());
                            if (isSameAnnotation(annotation.transcript(), featureProteinAnnotation, snpeffProteinAnnotation)) {
                                matchFound = true;
                            }
                        }
                    }

                    if (matchFound) {
                        LOGGER.debug("Found a match amongst candidate transcripts for '{}' on '{}", featureProteinAnnotation, featureGene);
                        matchCount++;
                    } else {
                        LOGGER.warn("Could not find a match amongst candidate transcripts for '{}' on '{}'",
                                featureProteinAnnotation,
                                featureGene);
                        diffCount++;
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

    @Nullable
    private static SnpEffAnnotation annotationForTranscript(@NotNull List<SnpEffAnnotation> annotations, @NotNull String transcript) {
        for (SnpEffAnnotation annotation : annotations) {
            if (annotation.isTranscriptFeature() && annotation.transcript().equals(transcript)) {
                return annotation;
            }
        }
        return null;
    }

    private boolean isSameAnnotation(@NotNull String transcript, @NotNull String featureAnnotation, @NotNull String snpeffAnnotation) {
        String curatedFeatureAnnotation = curateStartCodonAnnotation(featureAnnotation);
        if (curatedFeatureAnnotation.equals(snpeffAnnotation)) {
            return true;
        }

        Map<String, List<String>> transcriptMapping = SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT.get(transcript);
        if (transcriptMapping != null) {
            Set<String> requestedAnnotations = annotationsRequestedForMappingPerTranscript.get(transcript);
            if (requestedAnnotations == null) {
                requestedAnnotations = Sets.newHashSet(featureAnnotation);
            } else {
                requestedAnnotations.add(featureAnnotation);
            }
            annotationsRequestedForMappingPerTranscript.put(transcript, requestedAnnotations);
            List<String> mappedAnnotations = transcriptMapping.get(featureAnnotation);
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

        return serveToSnpEffMappings;
    }

    @NotNull
    private static Map<String, List<String>> createPDGFRAMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.D842_I843delinsIM", Lists.newArrayList("p.DI842IM"));
        map.put("p.D842_H845del", Lists.newArrayList("p.I843_D846del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createKITMap() {
        Map<String, List<String>> map = Maps.newHashMap();
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
        map.put("p.S768_V769insVAS", Lists.newArrayList("p.V769_D770insASV"));
        map.put("p.D770_N771insD", Lists.newArrayList("p.D770dup"));
        map.put("p.V769_D770insASV", Lists.newArrayList("p.A767_V769dup"));
        map.put("p.D770_N771insSVD", Lists.newArrayList("p.S768_D770dup"));
        map.put("p.D770_N771insNPG", Lists.newArrayList("p.P772_H773insGNP"));
        map.put("p.H773_V774insH", Lists.newArrayList("p.H773dup"));
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
        map.put("p.M774_A775insAYVM", Lists.newArrayList("p.A775_G776insYVMA"));
        map.put("p.A775_G776insYVMA", Lists.newArrayList("p.Y772_A775dup"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createBRCA1Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.V1688del", Lists.newArrayList("p.V1687del"));
        map.put("p.S1297del", Lists.newArrayList("p.S1298del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createARAFMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.Q347_A348del", Lists.newArrayList("p.Q349_A350del"));
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
}
