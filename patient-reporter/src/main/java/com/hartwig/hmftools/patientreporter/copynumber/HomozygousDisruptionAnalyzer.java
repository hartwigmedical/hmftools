package com.hartwig.hmftools.patientreporter.copynumber;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriverFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HomozygousDisruptionAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(HomozygousDisruptionAnalyzer.class);

    private HomozygousDisruptionAnalyzer() {

    }

    @NotNull
    public static List<LinxDriver> readingLinxDriver(@NotNull String linxDriversTsv) throws IOException {
        List<LinxDriver> linxDrivers = LinxDriverFile.read(linxDriversTsv);
        LOGGER.info("Loaded {} linx drivers from {}", linxDrivers.size(), linxDriversTsv);

        return extractDelDisruptions(linxDrivers);

    }

    @NotNull
    public static List<LinxDriver> extractDelDisruptions(@NotNull List<LinxDriver> linxDrivers) {
        List<LinxDriver> linxDriverList = Lists.newArrayList();
        for (LinxDriver homozygousDisruption : linxDrivers) {
            if (homozygousDisruption.eventType().equals("HOM_DEL_DISRUPTION")) {
                linxDriverList.add(ImmutableLinxDriver.builder()
                        .clusterId(homozygousDisruption.clusterId())
                        .gene(homozygousDisruption.gene())
                        .eventType(homozygousDisruption.eventType())
                        .build());
            }
        }
        return linxDriverList;
    }
}
