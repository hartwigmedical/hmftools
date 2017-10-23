package com.hartwig.hmftools.purple.segment;

import java.io.IOException;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.pcf.PCFFile;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RatioBreakPoints {
    private static final Logger LOGGER = LogManager.getLogger(RatioBreakPoints.class);

    @NotNull
    public static Multimap<String, GenomePosition> createRatioSegments(@NotNull final CommonConfig config) throws IOException {

        final String referenceFile = PCFFile.generateRatioFilename(config.cobaltDirectory(), config.refSample());
        LOGGER.info("Loading reference PCF segments from {}", referenceFile);
        final Multimap<String, GenomePosition> referenceBreakPoint = PCFFile.readPositions(config.windowSize(), referenceFile);

        final String tumorFile = PCFFile.generateRatioFilename(config.cobaltDirectory(), config.tumorSample());
        LOGGER.info("Loading tumor PCF segments from {}", tumorFile);
        final Multimap<String, GenomePosition> tumorBreakPoints = PCFFile.readPositions(config.windowSize(), tumorFile);

        LOGGER.info("Merging reference and tumor ratio break points");
        return GenomePositions.union(referenceBreakPoint, tumorBreakPoints);
    }
}
