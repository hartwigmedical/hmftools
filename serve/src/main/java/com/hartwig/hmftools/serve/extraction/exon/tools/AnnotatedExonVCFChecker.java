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

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE exon VCF checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        int totalCount = 0;
        int matchCount = 0;
        int diffCount = 0;

        String annotatedExonVcf = System.getProperty("user.home") + "/hmf/tmp/annotatedExons.vcf";

        LOGGER.info("Loading exons from '{}'", annotatedExonVcf);
        AbstractFeatureReader<VariantContext, LineIterator> reader =
                AbstractFeatureReader.getFeatureReader(annotatedExonVcf, new VCFCodec(), false);
        for (VariantContext variant : reader.iterator()) {
            totalCount++;

            String[] inputParts = variant.getAttributeAsString("input", Strings.EMPTY).split("\\|");
            String inputGene = inputParts[0];
            String inputTranscript = inputParts[1].equals("null") ? null : inputParts[1];
            int inputExonId = Integer.parseInt(inputParts[2]);

            List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(variant);
            if (isMatch(inputGene, inputTranscript, inputExonId, annotations)) {
                matchCount++;
            } else {
                diffCount++;
            }
        }

        LOGGER.info("Done comparing {} exons: {} matches and {} differences found.", totalCount, matchCount, diffCount);
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

    private static boolean isMatch(@NotNull String inputGene, @Nullable String inputTranscript, int inputExonId,
            @NotNull List<SnpEffAnnotation> annotations) {
        if (inputTranscript != null) {
            SnpEffAnnotation annotation = annotationForTranscript(annotations, inputTranscript);

            if (annotation != null) {
                int snpeffExonId = extractExonId(annotation.rank());
                if (inputExonId == snpeffExonId) {
                    LOGGER.debug("Identical on gene '{}': SERVE input exon id '{}' vs SnpEff exon id '{}'",
                            inputGene,
                            inputExonId,
                            snpeffExonId);
                    return true;
                } else {
                    LOGGER.warn("Difference on gene '{}': SERVE input exon id '{}' vs SnpEff exon id '{}'",
                            inputGene,
                            inputExonId,
                            snpeffExonId);
                    return false;
                }
            } else {
                LOGGER.warn("No suitable annotation found on gene '{}': SERVE input exon id '{}'", inputGene, inputExonId);
                return false;
            }
        } else {
            // In case input transcript is missing we try to match against any transcript.
            boolean matchFound = false;
            for (SnpEffAnnotation annotation : annotations) {
                if (annotation.isTranscriptFeature()) {
                    int snpeffExonId = extractExonId(annotation.rank());
                    if (inputExonId == snpeffExonId) {
                        matchFound = true;
                    }
                }
            }

            if (matchFound) {
                LOGGER.debug("Found a match amongst candidate transcripts for '{}' on '{}", inputExonId, inputGene);
                return true;
            } else {
                LOGGER.warn("Could not find a match amongst candidate transcripts '{}' for on '{}'", inputExonId, inputGene);
                return false;
            }
        }
    }

    private static int extractExonId(@NotNull String snpeffExonRank) {
        // Assume format is "{RANK}/{TOTAL}"
        return Integer.parseInt(snpeffExonRank.split("/")[0]);
    }
}
