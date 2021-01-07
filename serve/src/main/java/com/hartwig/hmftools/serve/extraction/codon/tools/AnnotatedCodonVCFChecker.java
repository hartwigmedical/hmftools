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

    private enum MatchType {
        IDENTICAL,
        NO_MATCH
    }

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE codon VCF checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        int totalCount = 0;
        int matchCount = 0;
        int diffCount = 0;

        String annotatedCodonVcf = System.getProperty("user.home") + "/hmf/tmp/annotated_codon.vcf";

        LOGGER.info("Loading codons from '{}'", annotatedCodonVcf);
        AbstractFeatureReader<VariantContext, LineIterator> reader =
                AbstractFeatureReader.getFeatureReader(annotatedCodonVcf, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            totalCount++;

            String[] inputParts = variant.getAttributeAsString("input", Strings.EMPTY).split("\\|");
            String inputGene = inputParts[0];
            String inputTranscript = inputParts[1].equals("null") ? null : inputParts[1];
            String inputProteinAnnotation = inputParts[2];
            String formattedHotspot = formatHotspot(variant);

            if (inputProteinAnnotation.equals(NO_INPUT_PROTEIN)) {
                LOGGER.debug("Skipping non-coding hotspot on '{}'", formattedHotspot);
            } else {
                List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(variant);
                MatchType match = determineMatch(inputGene, inputTranscript, inputProteinAnnotation, annotations);

                switch (match) {
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
        LOGGER.info("Done comparing {} records: {} matches and {} differences found.", totalCount, matchCount, diffCount);

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

    private static MatchType determineMatch(@NotNull String inputGene, @Nullable String inputTranscript,
            @NotNull String inputProteinAnnotation, @NotNull List<SnpEffAnnotation> annotations) {
        String inputCodon = inputProteinAnnotation.substring(0, inputProteinAnnotation.length() - 1);
        String snpEffCodon = Strings.EMPTY;
        String snpeffProteinAnnotation = Strings.EMPTY;

        if (inputTranscript != null) {
            SnpEffAnnotation annotation = annotationForTranscript(annotations, inputTranscript);

            if (annotation != null) {
                snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());
                snpEffCodon = snpeffProteinAnnotation.substring(0, snpeffProteinAnnotation.length() - 1);
                if (!isSameAnnotation(inputCodon, snpEffCodon)) {
                    LOGGER.warn("Difference on gene '{}' : SERVE input protein '{}' vs SnpEff protein '{}'",
                            inputGene,
                            inputCodon,
                            snpEffCodon);
                    return MatchType.NO_MATCH;
                } else {
                    LOGGER.debug("Identical on gene '{}' : SERVE input protein '{}' vs SnpEff protein '{}'",
                            inputGene,
                            inputCodon,
                            snpEffCodon);
                    return MatchType.IDENTICAL;
                }
            } else {
                LOGGER.warn("No match found on gene '{}': SERVE input protein '{}' vs SnpEff protein '{}'",
                        inputGene,
                        inputCodon,
                        snpEffCodon);
                return MatchType.NO_MATCH;
            }
        } else {
            // In case input transcript is missing we try to match against any transcript.
            boolean matchFound = false;
            for (SnpEffAnnotation annotation : annotations) {
                if (annotation.isTranscriptFeature()) {
                    snpeffProteinAnnotation = AminoAcidFunctions.forceSingleLetterProteinAnnotation(annotation.hgvsProtein());
                    snpEffCodon = snpeffProteinAnnotation.substring(0, snpeffProteinAnnotation.length() - 1);
                    if (isSameAnnotation(inputCodon, snpEffCodon)) {
                        matchFound = true;
                    }
                }
            }

            if (matchFound) {
                LOGGER.debug("Could not find a match amongst candidate transcripts for '{}' on '{}' of snpeff annotation '{}'",
                        inputProteinAnnotation,
                        inputGene,
                        snpeffProteinAnnotation);
                return MatchType.IDENTICAL;
            } else {
                LOGGER.warn("Found a match amongst candidate transcripts for '{}' on '{} of snpeff annotation '{}'",
                        inputProteinAnnotation,
                        inputGene,
                        snpeffProteinAnnotation);
                return MatchType.NO_MATCH;
            }
        }
    }

    private static boolean isSameAnnotation(@NotNull String inputAnnotation, @NotNull String snpeffAnnotation) {
        if (inputAnnotation.equals(snpeffAnnotation)) {
            return true;
        } else {
            return false;
        }
    }

    @NotNull
    private static String formatHotspot(@NotNull VariantContext variant) {
        return variant.getContig() + ":" + variant.getStart() + " " + variant.getReference().getBaseString() + ">"
                + variant.getAlternateAllele(0).getBaseString();
    }
}
