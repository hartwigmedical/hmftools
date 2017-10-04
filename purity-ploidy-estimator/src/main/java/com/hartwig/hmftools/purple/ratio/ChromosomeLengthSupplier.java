package com.hartwig.hmftools.purple.ratio;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengthFactory;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengthFile;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ChromosomeLengthSupplier implements Supplier<Map<String, ChromosomeLength>> {

    private static final Logger LOGGER = LogManager.getLogger(ChromosomeLengthSupplier.class);


    private final Map<String, ChromosomeLength> chromosomeLengths;

    public ChromosomeLengthSupplier(final CommonConfig config, Multimap<String, ReadRatio> tumorReadRatio) throws IOException, HartwigException {

        final String chrLengthFile = ChromosomeLengthFile.generateFilename(config.cobaltDirectory(), config.tumorSample());
        if (new File(chrLengthFile).exists()) {
            LOGGER.info("Loading chromosome lengths from {}", chrLengthFile);
            chromosomeLengths =
                    ChromosomeLengthFile.read(chrLengthFile).stream().collect(Collectors.toMap(ChromosomeLength::chromosome, item -> item));
        } else {
            LOGGER.info("Generating chromosome lengths from tumor read ratios");
            chromosomeLengths = ChromosomeLengthFactory.create(config.windowSize(), tumorReadRatio);
        }
    }

    @Override
    public Map<String, ChromosomeLength> get() {
        return chromosomeLengths;
    }
}
