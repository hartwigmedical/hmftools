package com.hartwig.hmftools.serve;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;
import com.hartwig.hmftools.serve.util.AminoAcidFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class ServeAnnotatedHotspotVCFChecker {

    private static final Logger LOGGER = LogManager.getLogger(ServeAnnotatedHotspotVCFChecker.class);

    private static final Map<String, Map<String, List<String>>> SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT = createMappings();

    public static void main(String[] args) throws IOException {
        //        Configurator.setRootLevel(Level.DEBUG);

        String annotatedHotspotVcf = System.getProperty("user.home") + "/hmf/tmp/annotatedHotspotsVicc.vcf";
        int totalCount = 0;
        int matchCount = 0;
        int whitelistedMatchCount = 0;
        int diffCount = 0;

        LOGGER.info("Loading hotspots from '{}'", annotatedHotspotVcf);
        AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(annotatedHotspotVcf, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            totalCount++;
            String[] featureParts = variant.getAttributeAsString("feature", Strings.EMPTY).split("\\|");
            String featureGene = featureParts[0];
            String featureTranscript = featureParts[1];
            String featureProteinAnnotation = featureParts[2];
            List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(variant);
            for (SnpEffAnnotation annotation : annotations) {
                if (annotation.isTranscriptFeature() && annotation.transcript().equals(featureTranscript)) {
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
                }
            }
        }

        LOGGER.info("Done comparing {} records: {} matches (of which {} through whitelisting) and {} differences found.",
                totalCount,
                matchCount,
                whitelistedMatchCount,
                diffCount);
    }

    private static boolean isSameAnnotation(@NotNull String transcript, @NotNull String featureAnnotation,
            @NotNull String snpeffAnnotation) {
        String curatedFeatureAnnotation = curateStartCodonAnnotation(featureAnnotation);
        if (curatedFeatureAnnotation.equals(snpeffAnnotation)) {
            return true;
        }

        Map<String, List<String>> transcriptMapping = SERVE_TO_SNPEFF_MAPPINGS_PER_TRANSCRIPT.get(transcript);
        if (transcriptMapping != null) {
            List<String> mappedAnnotations = transcriptMapping.get(featureAnnotation);
            if (mappedAnnotations != null) {
                return mappedAnnotations.contains(snpeffAnnotation);
            }
        }

        return false;
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
        serveToSnpEffMappings.put("ENST00000227507", createCCND1Map());
        serveToSnpEffMappings.put("ENST00000256078", createKRASMap());
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
        map.put("p.I559_D560insDKRMNS", Lists.newArrayList("p.Y556_K561dup"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createEGFRMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.S768_V769insVAS", Lists.newArrayList("p.A767_V769dup"));
        map.put("p.D770_N771insD", Lists.newArrayList("p.D770dup"));
        map.put("p.V769_D770insASV", Lists.newArrayList("p.A767_V769dup"));
        map.put("p.D770_N771insSVD", Lists.newArrayList("p.S768_D770dup"));
        map.put("p.D770_N771insNPG", Lists.newArrayList("p.D770_P772dup"));
        map.put("p.H773_V774insH", Lists.newArrayList("p.H773dup"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createBRAFMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.T599_V600insV", Lists.newArrayList("p.V600dup"));
        map.put("p.T599_V600insEAT", Lists.newArrayList("p.G596_A598dup"));
        map.put("p.T599_V600insETT", Lists.newArrayList("p.G596_A598dup"));
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
    private static Map<String, List<String>> createCCND1Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.L283_D294del", Lists.newArrayList("p.V281_D292del", "p.D282_V293del"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createKRASMap() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.A11_G12insGA", Lists.newArrayList("p.G10_A11dup"));
        return map;
    }

    @NotNull
    private static Map<String, List<String>> createFLT3Map() {
        Map<String, List<String>> map = Maps.newHashMap();
        map.put("p.S840_N841insGS", Lists.newArrayList("p.D839_S840dup"));
        map.put("p.K602_W603insYEYDLK", Lists.newArrayList("p.Y597_K602dup"));
        map.put("p.W603_E604insDREYEYDLKW", Lists.newArrayList("p.D593_K602dup"));
        map.put("p.L601_K602insREYEYDL", Lists.newArrayList("p.R595_L601dup"));
        map.put("p.D600_L601insFREYEYD", Lists.newArrayList("p.F594_D600dup"));
        map.put("p.Y599_D600insPAPQIMSTSTLISENMNIA", Lists.newArrayList("p.V581_Y599dup"));
        map.put("p.Y599_D600insEYEYEYEY", Lists.newArrayList("p.Y591_E598dup"));
        map.put("p.Y599_D600insGLYVDFREYEY", Lists.newArrayList("p.E588_E598dup"));
        map.put("p.D600_L601insDFREYEYD", Lists.newArrayList("p.D593_D600dup"));
        map.put("p.E598_Y599insDVDFREYE", Lists.newArrayList("p.Y591_E598dup"));
        map.put("p.Y599_D600insSTDNEYFYVDFREYEY", Lists.newArrayList("p.G583_E598dup"));
        map.put("p.E598_Y599insGLVQVTGSSDNEYFYVDFREYE", Lists.newArrayList("p.Q577_E598dup"));
        map.put("p.F594_R595insSDNEYFYVDF", Lists.newArrayList("p.S585_F594dup"));
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
        map.put("p.M774_A775insAYVM", Lists.newArrayList("p.Y772_A775dup"));
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
