package com.hartwig.hmftools.serve.extraction.exon.tools;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;
import com.hartwig.hmftools.serve.extraction.codon.tools.AnnotatedCodonVCFChecker;

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

public class AnnotatedExonVCFChecker {

    private static final Logger LOGGER = LogManager.getLogger(AnnotatedCodonVCFChecker.class);
    private static final boolean LOG_DEBUG = true;

    private static final String NO_INPUT_PROTEIN = "-";

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE codon VCF checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String annotatedCodonVcf = System.getProperty("user.home") + "/hmf/tmp/annotated.vcf";

        AbstractFeatureReader<VariantContext, LineIterator> reader =
                AbstractFeatureReader.getFeatureReader(annotatedCodonVcf, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            String[] inputParts = variant.getAttributeAsString("input", Strings.EMPTY).split("\\|");
            String inputGene = inputParts[0];
            String inputTranscript = inputParts[1].equals("null") ? null : inputParts[1];
            String inputExonId = inputParts[2];

            List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(variant);
            determineMatch(inputGene, inputTranscript, inputExonId, annotations);

        }
    }

    @Nullable
    private static SnpEffAnnotation annotationForTranscript(@NotNull List<SnpEffAnnotation> annotations, @NotNull String transcript) {
        for (SnpEffAnnotation annotation : annotations) {
            if (annotation.isTranscriptFeature()) {
                return annotation;
            }
        }
        return null;
    }

    private static void determineMatch(@NotNull String inputGene, @Nullable String inputTranscript, @NotNull String inputExonId,
            @NotNull List<SnpEffAnnotation> annotations) {
        String snpeffExonID = Strings.EMPTY;
        if (inputTranscript != null) {
            SnpEffAnnotation annotation = annotationForTranscript(annotations, inputTranscript);

            if (annotation != null) {
                snpeffExonID = annotation.rank();

            } else {
                LOGGER.warn("Could not find snpeff annotation for '{}' on '{}'!", inputTranscript, inputGene);
            }
        } else {
            // In case input transcript is missing we try to match against any transcript.
            for (SnpEffAnnotation annotation : annotations) {
                if (annotation.isTranscriptFeature()) {
                    snpeffExonID = annotation.rank();

                }
            }
        }
        String snpeffExonIDExtract = snpeffExonID.split("/")[0];
        if (inputExonId.equals(snpeffExonIDExtract)) {
            LOGGER.debug("Found a match amongst candidate transcripts for '{}' on '{} of snpeff annotation '{}'",
                    inputExonId,
                    inputGene,
                    snpeffExonIDExtract);
        } else {
            LOGGER.warn("Could not find a match amongst candidate transcripts '{}' for on '{}' of snpeff annotation '{}'",
                    inputExonId,
                    inputGene,
                    snpeffExonIDExtract);
        }
    }

}
