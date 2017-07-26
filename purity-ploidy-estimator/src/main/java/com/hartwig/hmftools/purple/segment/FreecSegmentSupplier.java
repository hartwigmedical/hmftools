package com.hartwig.hmftools.purple.segment;

import java.util.List;
import java.util.function.Supplier;

import com.hartwig.hmftools.common.copynumber.freec.FreecRatioRegions;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.purple.ratio.FreecRatioSupplier;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class FreecSegmentSupplier implements Supplier<List<GenomeRegion>> {

    private static final Logger LOGGER = LogManager.getLogger(FreecSegmentSupplier.class);

    private final List<GenomeRegion> segments;

    public FreecSegmentSupplier(@NotNull final FreecRatioSupplier suppler) {
        LOGGER.info("Using freec segmentation");
        segments = FreecRatioRegions.createRegionsFromRatios(suppler.tumorFreecRatios());
    }

    @Override
    public List<GenomeRegion> get() {
        return segments;
    }
}
