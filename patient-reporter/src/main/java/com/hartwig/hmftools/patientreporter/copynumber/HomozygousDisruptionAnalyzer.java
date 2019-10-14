package com.hartwig.hmftools.patientreporter.copynumber;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.HomozygousDisruptionFile;
import com.hartwig.hmftools.common.linx.ImmutableHomozygousDisruption;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HomozygousDisruptionAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(HomozygousDisruptionAnalyzer.class);

    private HomozygousDisruptionAnalyzer() {

    }

    @NotNull
    public static List<HomozygousDisruption> analysisHomozygousDisruption(@NotNull String linxDriversTsv) throws IOException {
        List<HomozygousDisruption> linxDrivers = HomozygousDisruptionFile.read(linxDriversTsv);
        LOGGER.info("Loaded {} linx drivers from {}", linxDrivers.size(), linxDriversTsv);

        return extractDelDisruptions(linxDrivers);

    }

    @NotNull
    public static List<HomozygousDisruption> extractDelDisruptions(@NotNull List<HomozygousDisruption> linxDrivers) {
        List<HomozygousDisruption> homozygousDisruptionList = Lists.newArrayList();
        for (HomozygousDisruption homozygousDisruption : linxDrivers) {
            if (homozygousDisruption.eventType().equals("HOM_DEL_DISRUPTION")) {
                homozygousDisruptionList.add(ImmutableHomozygousDisruption.builder()
                        .clusterId(homozygousDisruption.clusterId())
                        .gene(homozygousDisruption.gene())
                        .eventType(homozygousDisruption.eventType())
                        .build());
            }
        }
        return homozygousDisruptionList;
    }
}
