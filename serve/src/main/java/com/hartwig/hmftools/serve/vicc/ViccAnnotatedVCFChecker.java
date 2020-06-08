package com.hartwig.hmftools.serve.vicc;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;
import com.hartwig.hmftools.serve.util.AminoAcidFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class ViccAnnotatedVCFChecker {

    private static final Logger LOGGER = LogManager.getLogger(ViccAnnotatedVCFChecker.class);

    public static void main(String[] args) throws IOException {
        String annotatedHotspotVcf = System.getProperty("user.home") + "/hmf/tmp/annotatedHotspotsVicc.vcf";

        AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(annotatedHotspotVcf, new VCFCodec(), false);

        for (VariantContext variant : reader.iterator()) {
            String[] featureParts = variant.getAttributeAsString("feature", Strings.EMPTY).split("\\|");
            String featureGene = featureParts[0];
            String featureTranscript = featureParts[1];
            String featureProteinAnnotation = featureParts[2];

            List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(variant);
            for (SnpEffAnnotation annotation : annotations) {
                if (annotation.isTranscriptFeature() && annotation.transcript().equals(featureTranscript)) {
                    String snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());
                    if (!snpeffProteinAnnotation.equals(featureProteinAnnotation)) {
                        LOGGER.info("Difference on gene '{}' - {}:{} {}->{} : Feature protein '{}' vs SnpEff protein '{}'",
                                featureGene,
                                variant.getContig(),
                                variant.getStart(),
                                variant.getReference().getBaseString(),
                                variant.getAlternateAllele(0).getBaseString(),
                                featureProteinAnnotation,
                                snpeffProteinAnnotation);
                    }
                }
            }
        }
    }
}
