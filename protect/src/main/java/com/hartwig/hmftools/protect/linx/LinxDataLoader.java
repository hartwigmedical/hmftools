package com.hartwig.hmftools.protect.linx;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
import com.hartwig.hmftools.protect.ProtectConfig;
import com.hartwig.hmftools.protect.homozygousdisruption.HomozygousDisruptionAnalyzer;
import com.hartwig.hmftools.protect.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.protect.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.protect.structural.ReportableGeneDisruptionFactory;
import com.hartwig.hmftools.protect.viralinsertion.ViralInsertion;
import com.hartwig.hmftools.protect.viralinsertion.ViralInsertionAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LinxDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(PurpleDataLoader.class);

    @NotNull
    public static LinxData load(ProtectConfig config) throws IOException {
        return load(config.linxFusionTsv(), config.linxBreakendTsv(), config.linxViralInsertionTsv(), config.linxDriversTsv());
    }

    @NotNull
    public static LinxData load(@NotNull String linxFusionTsv, @NotNull String linxBreakendTsv, @NotNull String linxViralInsertionTsv,
            @NotNull String linxDriversTsv) throws IOException {
        LOGGER.info("Loaded LINX data from {}", new File(linxFusionTsv).getParent());
        List<LinxFusion> linxFusions = LinxFusion.read(linxFusionTsv).stream().filter(LinxFusion::reported).collect(Collectors.toList());
        LOGGER.info(" Reportable fusions: {}", linxFusions.size());

        List<LinxBreakend> linxBreakends =
                LinxBreakend.read(linxBreakendTsv).stream().filter(LinxBreakend::reportedDisruption).collect(Collectors.toList());
        List<ReportableGeneDisruption> reportableGeneDisruptions = ReportableGeneDisruptionFactory.convert(linxBreakends);
        LOGGER.info(" Reportable disruptions: {}", reportableGeneDisruptions.size());

        List<LinxViralInsertion> viralInsertionList = LinxViralInsertion.read(linxViralInsertionTsv);
        List<ViralInsertion> reportableViralInsertions = ViralInsertionAnalyzer.analyzeViralInsertions(viralInsertionList);
        LOGGER.info(" Reportable viral insertions: {}", reportableViralInsertions.size());

        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions =
                HomozygousDisruptionAnalyzer.extractFromLinxDriversTsv(linxDriversTsv);
        LOGGER.info(" Reportable homozygous disruptions: {}", reportableHomozygousDisruptions.size());

        return ImmutableLinxData.builder()
                .addAllFusions(linxFusions)
                .addAllGeneDisruptions(reportableGeneDisruptions)
                .addAllViralInsertions(reportableViralInsertions)
                .addAllHomozygousDisruptions(reportableHomozygousDisruptions)
                .build();
    }
}
