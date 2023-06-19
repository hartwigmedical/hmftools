package com.hartwig.hmftools.datamodel.orange;

import java.time.LocalDate;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeConfig {

    Logger LOGGER = LogManager.getLogger(OrangeConfig.class);

    @NotNull
    String tumorSampleId();

    @Nullable
    String referenceSampleId();

    @Nullable
    OrangeRNAConfig rnaConfig();

    @NotNull
    Set<String> primaryTumorDoids();

    @NotNull
    LocalDate experimentDate();

    @NotNull
    OrangeRefGenomeVersion refGenomeVersion();

    @NotNull
    String outputDir();

    @NotNull
    String doidJsonFile();

    @NotNull
    String cohortMappingTsv();

    @NotNull
    String cohortPercentilesTsv();

    @NotNull
    String driverGenePanelTsv();

    @NotNull
    String knownFusionFile();

    @NotNull
    String ensemblDataDirectory();

    @Nullable
    String pipelineVersionFile();

    @Nullable
    String refSampleWGSMetricsFile();

    @Nullable
    String refSampleFlagstatFile();

    @NotNull
    String tumorSampleWGSMetricsFile();

    @NotNull
    String tumorSampleFlagstatFile();

    @Nullable
    String sageGermlineGeneCoverageTsv();

    @Nullable
    String sageSomaticRefSampleBQRPlot();

    @NotNull
    String sageSomaticTumorSampleBQRPlot();

    @NotNull
    String purpleDataDirectory();

    @NotNull
    String purplePlotDirectory();

    @NotNull
    String linxSomaticDataDirectory();

    @Nullable
    String linxGermlineDataDirectory();

    @NotNull
    String linxPlotDirectory();

    @NotNull
    String lilacResultCsv();

    @NotNull
    String lilacQcCsv();

    @Nullable
    String annotatedVirusTsv();

    @Nullable
    String chordPredictionTxt();

    @Nullable
    String cuppaResultCsv();

    @Nullable
    String cuppaSummaryPlot();

    @Nullable
    String cuppaFeaturePlot();

    @Nullable
    String cuppaChartPlot();

    @Nullable
    String peachGenotypeTsv();

    @Nullable
    String sigsAllocationTsv();

    boolean convertGermlineToSomatic();

    boolean limitJsonOutput();

    boolean addDisclaimer();

}
