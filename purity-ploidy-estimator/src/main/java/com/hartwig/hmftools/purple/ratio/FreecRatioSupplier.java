package com.hartwig.hmftools.purple.ratio;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatio;
import com.hartwig.hmftools.common.copynumber.freec.FreecRatioFactory;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.ratio.ImmutableReadRatio;
import com.hartwig.hmftools.common.purple.ratio.ReadRatio;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class FreecRatioSupplier implements RatioSupplier {

    private static final Logger LOGGER = LogManager.getLogger(FreecRatioSupplier.class);

    private final Multimap<String, ReadRatio> referenceRatios;
    private final Multimap<String, ReadRatio> tumorRatios;
    private final List<FreecRatio> tumorFreecRatios;

    public FreecRatioSupplier(final CommonConfig config) throws IOException, HartwigException {

        // KODU: Even though this retrieves normal ratios, freec uses the tumor sample name in the file name.
        LOGGER.info("Reading freec reference ratios");
        referenceRatios = toRatio(FreecRatioFactory.loadNormalRatios(config.freecDirectory(), config.tumorSample()));

        LOGGER.info("Reading freec tumor ratios");
        tumorFreecRatios = FreecRatioFactory.loadTumorRatios(config.freecDirectory(), config.tumorSample());
        tumorRatios = toRatio(tumorFreecRatios);
    }

    @Override
    @NotNull
    public Multimap<String, ReadRatio> referenceRatios() {
        return referenceRatios;
    }

    @Override
    @NotNull
    public Multimap<String, ReadRatio> tumorRatios() {
        return tumorRatios;
    }

    public List<FreecRatio> tumorFreecRatios() {
        return tumorFreecRatios;
    }

    private Multimap<String, ReadRatio> toRatio(List<FreecRatio> ratios) {
        final Multimap<String, ReadRatio> results = ArrayListMultimap.create();
        for (FreecRatio ratio : ratios) {
            results.put(ratio.chromosome(), ImmutableReadRatio.builder().from(ratio).ratio(ratio.ratio()).build());
        }

        return results;
    }
}
