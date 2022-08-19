package com.hartwig.hmftools.serve.extraction.codon.tools;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.serve.extraction.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.serve.extraction.snpeff.SnpEffAnnotationParser;
import com.hartwig.hmftools.serve.extraction.util.VCFWriterFactory;

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

    private static final boolean LOG_DEBUG = false;

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE Codon VCF checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        int totalCount = 0;
        int matchCount = 0;
        int diffCount = 0;

        String annotatedCodonVcf = System.getProperty("user.home") + "/hmf/tmp/annotatedCodons.vcf";

        LOGGER.info("Loading codons from '{}'", annotatedCodonVcf);
        AbstractFeatureReader<VariantContext, LineIterator> reader =
                AbstractFeatureReader.getFeatureReader(annotatedCodonVcf, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            totalCount++;

            String[] inputParts = variant.getAttributeAsString(VCFWriterFactory.INPUT_FIELD, Strings.EMPTY).split("\\|");
            String inputGene = inputParts[0];
            String inputTranscript = inputParts[1].equals("null") ? null : inputParts[1];
            int inputCodon = Integer.parseInt(inputParts[2]);

            List<SnpEffAnnotation> annotations = SnpEffAnnotationParser.fromContext(variant);
            if (isMatch(inputGene, inputTranscript, inputCodon, annotations)) {
                matchCount++;
            } else {
                diffCount++;
            }
        }

        LOGGER.info("Done comparing {} codons: {} matches and {} differences found.", totalCount, matchCount, diffCount);
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

    private static boolean isMatch(@NotNull String inputGene, @Nullable String inputTranscript, int inputCodon,
            @NotNull List<SnpEffAnnotation> annotations) {
        if (inputTranscript != null) {
            SnpEffAnnotation annotation = annotationForTranscript(annotations, inputTranscript);

            if (annotation != null) {
                int snpEffCodon = extractCodon(annotation.hgvsProtein());
                if (inputCodon == snpEffCodon) {
                    LOGGER.debug("Identical on gene '{}': SERVE input codon '{}' vs SnpEff codon '{}'", inputGene, inputCodon, snpEffCodon);
                    return true;
                } else {
                    LOGGER.warn("Difference on gene '{}': SERVE input codon '{}' vs SnpEff codon '{}'", inputGene, inputCodon, snpEffCodon);
                    return false;
                }
            } else {
                LOGGER.warn("No suitable SnpEff annotation found on gene '{}': SERVE input codon '{}'", inputGene, inputCodon);
                return false;
            }
        } else {
            // In case input transcript is missing we try to match against any transcript.
            boolean matchFound = false;
            for (SnpEffAnnotation annotation : annotations) {
                if (annotation.isTranscriptFeature()) {
                    int snpEffCodon = extractCodon(annotation.hgvsProtein());
                    if (inputCodon == snpEffCodon) {
                        matchFound = true;
                    }
                }
            }

            if (matchFound) {
                LOGGER.debug("Found a match amongst candidate transcripts for '{}' on '{}", inputCodon, inputGene);
                return true;
            } else {
                LOGGER.warn("Could not find a match amongst candidate transcripts for '{}' on '{}'", inputCodon, inputGene);
                return false;
            }
        }
    }

    private static int extractCodon(@NotNull String hgvsProteinAnnotation) {
        String singleLetterAA = AminoAcids.forceSingleLetterProteinAnnotation(hgvsProteinAnnotation);
        // The single letter AA should always start with "p.{A}"
        return Integer.parseInt(singleLetterAA.substring(3, singleLetterAA.length() - 1));
    }
}
