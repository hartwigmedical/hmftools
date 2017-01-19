package com.hartwig.hmftools.patientreporter.slicing;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HMFSlicingAnnotation {

    private static final Logger LOGGER = LogManager.getLogger(HMFSlicingAnnotation.class);

    // KODU: We expect a format like "ENST00001.1 (GEN)"
    private static final String TRANSCRIPT_GENE_SEPARATOR = " ";
    private static final String TRANSCRIPT_VERSION_SEPARATOR = "\\.";

    @NotNull
    private final String transcriptID;
    private final int transcriptVersion;
    @NotNull
    private final String gene;

    @Nullable
    public static HMFSlicingAnnotation fromGenomeRegion(@NotNull final GenomeRegion region) {
        final String annotation = region.annotation();
        if (annotation != null) {
            final String[] parts = annotation.split(TRANSCRIPT_GENE_SEPARATOR);
            if (parts.length == 2 && parts[1].trim().length() > 2) {
                final String[] transcriptParts = parts[0].trim().split(TRANSCRIPT_VERSION_SEPARATOR);
                if (transcriptParts.length == 2) {
                    final String gene = parts[1].trim().substring(1, parts[1].length() - 1);
                    return new HMFSlicingAnnotation(transcriptParts[0], Integer.valueOf(transcriptParts[1]), gene);
                } else {
                    LOGGER.warn("Transcript part does not have correct format: " + parts[0].trim());
                }
            } else {
                LOGGER.warn("Annotation part does not have correct format: " + annotation);
            }
        } else {
            LOGGER.warn("No annotation present on " + region);
        }

        return null;
    }

    private HMFSlicingAnnotation(@NotNull final String transcriptID, final int transcriptVersion,
            @NotNull final String gene) {
        this.transcriptID = transcriptID;
        this.transcriptVersion = transcriptVersion;
        this.gene = gene;
    }

    @NotNull
    public String transcriptID() {
        return transcriptID;
    }

    int transcriptVersion() {
        return transcriptVersion;
    }

    @NotNull
    String gene() {
        return gene;
    }
}
