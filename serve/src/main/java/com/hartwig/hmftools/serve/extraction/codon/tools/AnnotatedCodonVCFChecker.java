package com.hartwig.hmftools.serve.extraction.codon.tools;

import java.io.IOException;
import java.util.List;

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

public class AnnotatedCodonVCFChecker {
    private static final Logger LOGGER = LogManager.getLogger(AnnotatedCodonVCFChecker.class);
    private static final boolean LOG_DEBUG = true;

    private static final String NO_INPUT_PROTEIN = "-";

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE codon VCF checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String annotatedCodonVcf = System.getProperty("user.home") + "/hmf/tmp/annotated_codon.vcf";

        AbstractFeatureReader<VariantContext, LineIterator> reader =
                AbstractFeatureReader.getFeatureReader(annotatedCodonVcf, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            String[] inputParts = variant.getAttributeAsString("input", Strings.EMPTY).split("\\|");
            String inputGene = inputParts[0];
            String inputTranscript = inputParts[1].equals("null") ? null : inputParts[1];
            String inputProteinAnnotation = inputParts[2];
            String formattedHotspot = formatHotspot(variant);

            if (inputProteinAnnotation.equals(NO_INPUT_PROTEIN)) {
                LOGGER.debug("Skipping non-coding hotspot on '{}'", formattedHotspot);
            } else {
                List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(variant);
                determineMatch(inputGene, inputTranscript, inputProteinAnnotation, annotations);
            }
        }
        LOGGER.info("Codons are checked!");
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

    private static void determineMatch(@NotNull String inputGene, @Nullable String inputTranscript, @NotNull String inputProteinAnnotation,
            @NotNull List<SnpEffAnnotation> annotations) {
        String snpeffProteinAnnotation = Strings.EMPTY;
        if (inputTranscript != null) {
            SnpEffAnnotation annotation = annotationForTranscript(annotations, inputTranscript);

            if (annotation != null) {
                snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());

            } else {
                LOGGER.warn("Could not find snpeff annotation for '{}' on '{}'!", inputTranscript, inputGene);
            }
        } else {
            // In case input transcript is missing we try to match against any transcript.
            for (SnpEffAnnotation annotation : annotations) {
                if (annotation.isTranscriptFeature()) {
                    snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());

                }
            }
        }
        String inputCodon = inputProteinAnnotation.substring(0, inputProteinAnnotation.length() - 1);
        String snpEffCodon = snpeffProteinAnnotation.substring(0, snpeffProteinAnnotation.length() - 1);

        if (inputCodon.equals(snpEffCodon)) {
            LOGGER.debug("Found a match amongst candidate transcripts for '{}' on '{} of snpeff annotation '{}'",
                    inputProteinAnnotation,
                    inputGene,
                    snpeffProteinAnnotation);
        } else {
            LOGGER.warn("Could not find a match amongst candidate transcripts for '{}' on '{}' of snpeff annotation '{}'",
                    inputProteinAnnotation,
                    inputGene,
                    snpeffProteinAnnotation);
        }
    }

    @NotNull
    private static String formatHotspot(@NotNull VariantContext variant) {
        return variant.getContig() + ":" + variant.getStart() + " " + variant.getReference().getBaseString() + ">"
                + variant.getAlternateAllele(0).getBaseString();
    }
}
