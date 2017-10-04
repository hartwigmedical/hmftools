package com.hartwig.hmftools.purple.ratio;

import java.io.IOException;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatioFile;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ReadCountRatioSupplier implements RatioSupplier {

    private static final Logger LOGGER = LogManager.getLogger(ReadCountRatioSupplier.class);

    private final Multimap<String, ReadRatio> referenceRatios;
    private final Multimap<String, ReadRatio> tumorRatios;

    public ReadCountRatioSupplier(final CommonConfig config)
            throws IOException, HartwigException {

        final String tumorRatioFile = ReadRatioFile.generateFilename(config.cobaltDirectory(), config.tumorSample());
        final String referenceRatioFile = ReadRatioFile.generateFilename(config.cobaltDirectory(), config.refSample());

        LOGGER.info("Loading reference ratios from {}", referenceRatioFile);
        referenceRatios = ReadRatioFile.read(referenceRatioFile);

        LOGGER.info("Loading tumor ratios from {}", tumorRatioFile);
        tumorRatios = ReadRatioFile.read(tumorRatioFile);
    }

    @Override
    @NotNull
    public Multimap<String, ReadRatio> tumorRatios() {
        return tumorRatios;
    }

    @Override
    @NotNull
    public Multimap<String, ReadRatio> referenceRatios() {
        return referenceRatios;
    }

}
