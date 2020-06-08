package com.hartwig.hmftools.serve.vicc;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;
import com.hartwig.hmftools.serve.util.AminoAcidFunctions;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class ViccAnnotatedVCFChecker {

    public static void main(String[] args) throws IOException {
        String annotatedHotspotVcf = System.getProperty("user.home") + "/hmf/tmp/annotatedHotspotsVicc.vcf";

        AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(annotatedHotspotVcf, new VCFCodec(), false);

        for (VariantContext variant : reader.iterator()) {
            String feature = variant.getAttributeAsString("feature", Strings.EMPTY);
            String featureGene = extractGene(feature);
            String featureTranscript = extractTranscript(feature);
            String featureProteinAnnotation = extractProteinAnnotation(feature);

            List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(variant);
            for (SnpEffAnnotation annotation : annotations) {
                if (annotation.isTranscriptFeature() && annotation.transcript().equals(featureTranscript)) {
                    String snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());
                }
            }
        }
    }

    @NotNull
    static String extractGene(@NotNull String feature) {
        return feature.substring(0, feature.indexOf(":"));
    }

    @NotNull
    static String extractProteinAnnotation(@NotNull String feature) {
        return feature.substring(feature.indexOf(":") + 1, feature.indexOf(" - "));
    }

    @NotNull
    static String extractTranscript(@NotNull String feature) {
        return feature.substring(feature.indexOf(" - ") + 3);
    }
}
