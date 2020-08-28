package com.hartwig.hmftools.serve.hotspot;

import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_DELETION;
import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_DUPLICATION;
import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_FRAMESHIFT_SUFFIX;
import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_FRAMESHIFT_SUFFIX_WITH_STOP_GAINED;
import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_INSERTION;
import static com.hartwig.hmftools.serve.util.HgvsConstants.HGVS_RANGE_INDICATOR;

import java.io.FileNotFoundException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.transvar.Transvar;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HotspotGenerator {

    private static final Logger LOGGER = LogManager.getLogger(HotspotGenerator.class);

    private static final int MAX_INFRAME_BASE_LENGTH = 50;

    @NotNull
    private final ProteinResolver proteinResolver;

    @NotNull
    public static HotspotGenerator transvarWithRefGenome(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile)
            throws FileNotFoundException {
        LOGGER.info("Creating hotspot generator with ref genome version '{}' and fasta path '{}'", refGenomeVersion, refGenomeFastaFile);
        return new HotspotGenerator(Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile));
    }

    @NotNull
    public static HotspotGenerator dummy() {
        return new HotspotGenerator(new ProteinResolver() {
            @NotNull
            @Override
            public List<VariantHotspot> extractHotspotsFromProteinAnnotation(@NotNull final String gene,
                    @Nullable final String specificTranscript, @NotNull final String proteinAnnotation) {
                return Lists.newArrayList();
            }

            @NotNull
            @Override
            public Set<String> unresolvedProteinAnnotations() {
                return Sets.newHashSet();
            }
        });
    }

    private HotspotGenerator(@NotNull ProteinResolver proteinResolver) {
        this.proteinResolver = proteinResolver;
    }

    @NotNull
    public List<VariantHotspot> generateHotspots(@NotNull String gene, @Nullable String specificTranscript,
            @NotNull String proteinAnnotation) {
        if (isResolvableProteinAnnotation(proteinAnnotation)) {
            return proteinResolver.extractHotspotsFromProteinAnnotation(gene,
                    specificTranscript,
                    proteinAnnotation);
        }

        return Lists.newArrayList();
    }

    @NotNull
    public Set<String> unresolvedProteinAnnotations() {
        return proteinResolver.unresolvedProteinAnnotations();
    }

    public static boolean isResolvableProteinAnnotation(@NotNull String proteinAnnotation) {
        try {
            if (isFrameshift(proteinAnnotation)) {
                return isValidFrameshift(proteinAnnotation);
            } else if (proteinAnnotation.contains(HGVS_RANGE_INDICATOR)) {
                return isValidRangeMutation(proteinAnnotation);
            } else if (proteinAnnotation.contains(HGVS_DELETION + HGVS_INSERTION)) {
                return isValidComplexDeletionInsertion(proteinAnnotation);
            } else if (proteinAnnotation.startsWith("*")) {
                return true;
            } else {
                return isValidSingleCodonMutation(proteinAnnotation);
            }
        } catch (Exception exception) {
            LOGGER.warn("Could not determine whether protein annotation is resolvable due to '{}'", exception.getMessage(), exception);
            return false;
        }
    }

    private static boolean isFrameshift(@NotNull String proteinAnnotation) {
        return proteinAnnotation.endsWith(HGVS_FRAMESHIFT_SUFFIX) || proteinAnnotation.endsWith(HGVS_FRAMESHIFT_SUFFIX_WITH_STOP_GAINED);
    }

    private static boolean isValidFrameshift(@NotNull String proteinAnnotation) {
        int frameshiftPosition = proteinAnnotation.indexOf(HGVS_FRAMESHIFT_SUFFIX);
        if (frameshiftPosition > 1) {
            return isInteger(proteinAnnotation.substring(frameshiftPosition - 1, frameshiftPosition));
        }

        return false;
    }

    private static boolean isValidRangeMutation(@NotNull String proteinAnnotation) {
        assert proteinAnnotation.contains(HGVS_RANGE_INDICATOR);

        // Features could be ranges such as E102_I103del. We whitelist specific feature types when analyzing a range.
        String[] annotationParts = proteinAnnotation.split(HGVS_RANGE_INDICATOR);
        String annotationStartPart = annotationParts[0];
        String annotationEndPart = annotationParts[1];
        if (annotationEndPart.contains(HGVS_INSERTION) || annotationEndPart.contains(HGVS_DUPLICATION)
                || annotationEndPart.contains(HGVS_DELETION)) {
            int indexOfEvent;
            // Keep in mind that 'del' always comes prior to 'ins' in situations of complex inframes.
            if (annotationEndPart.contains(HGVS_DELETION)) {
                indexOfEvent = annotationEndPart.indexOf(HGVS_DELETION);
            } else if (annotationEndPart.contains(HGVS_DUPLICATION)) {
                indexOfEvent = annotationEndPart.indexOf(HGVS_DUPLICATION);
            } else {
                indexOfEvent = annotationEndPart.indexOf(HGVS_INSERTION);
            }

            long start = Long.parseLong(annotationStartPart.substring(1));
            long end = Long.parseLong(annotationEndPart.substring(1, indexOfEvent));
            return 3 * (1 + end - start) <= MAX_INFRAME_BASE_LENGTH;
        } else {
            return false;
        }
    }

    private static boolean isValidComplexDeletionInsertion(@NotNull String proteinAnnotation) {
        String[] annotationParts = proteinAnnotation.split(HGVS_DELETION + HGVS_INSERTION);

        return isInteger(annotationParts[0].substring(1)) && (3 * annotationParts[1].length()) <= MAX_INFRAME_BASE_LENGTH;
    }

    private static boolean isValidSingleCodonMutation(@NotNull String proteinAnnotation) {
        if (proteinAnnotation.contains(HGVS_INSERTION)) {
            // Insertions are only allowed in a range, since we need to know where to insert the sequence exactly.
            return false;
        }

        // Features are expected to look something like V600E (1 char - N digits - M chars)
        if (proteinAnnotation.length() < 3) {
            return false;
        }

        if (!Character.isLetter(proteinAnnotation.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(proteinAnnotation.charAt(1))) {
            return false;
        }

        boolean haveObservedNonDigit = !Character.isDigit(proteinAnnotation.charAt(2));
        int firstNotDigit = haveObservedNonDigit ? 2 : -1;
        for (int i = 3; i < proteinAnnotation.length(); i++) {
            char charToEvaluate = proteinAnnotation.charAt(i);
            if (haveObservedNonDigit && Character.isDigit(charToEvaluate)) {
                return false;
            }
            boolean isDigit = Character.isDigit(charToEvaluate);
            if (!isDigit && firstNotDigit == -1) {
                firstNotDigit = i;
            }

            haveObservedNonDigit = haveObservedNonDigit || !isDigit;
        }

        if (!haveObservedNonDigit) {
            return false;
        }

        String newAminoAcid = proteinAnnotation.substring(firstNotDigit);
        // X is a wildcard that we don't support, and "/" indicates logical OR that we don't support.
        return !newAminoAcid.equals("X") && !newAminoAcid.contains("/");
    }

    private static boolean isInteger(@NotNull String value) {
        try {
            Integer.parseInt(value);
            return true;
        } catch (NumberFormatException exp) {
            return false;
        }
    }
}
