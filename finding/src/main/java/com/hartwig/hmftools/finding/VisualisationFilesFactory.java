package com.hartwig.hmftools.finding;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.datamodel.finding.VisualisationFiles;
import com.hartwig.hmftools.datamodel.finding.VisualisationFilesBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jspecify.annotations.Nullable;

public class VisualisationFilesFactory
{
    private static final Logger LOGGER = LogManager.getLogger(VisualisationFilesFactory.class);
    private final String visDir;

    public VisualisationFilesFactory(String visDir)
    {
        this.visDir = visDir;
        LOGGER.info("Creating visualisation directory: {}", visDir);
        FileWriterUtils.checkCreateOutputDir(visDir);
    }

    public VisualisationFiles create(
            final String pipelineOutputDir,
            final String tumorSampleId,
            @Nullable final String referenceSampleId,
            boolean hasCuppa)
    {
        LOGGER.info("Loading plots");

        List<String> linxDriverPlots = new ArrayList<>();
        Path linxPlotDir = Paths.get(pipelineOutputDir + File.separator + "linx/somatic_plots/reportable/");
        if (Files.exists(linxPlotDir) && Files.isDirectory(linxPlotDir))
        {
            try
            {
                try (Stream<Path> stream = Files.list(linxPlotDir))
                {
                    stream.forEach(p -> {
                        linxDriverPlots.add(processVisualisation(p.toString()));
                    });
                }
            } catch (IOException e) {
                throw new RuntimeException("Failed to list LINX plot directory: " + linxPlotDir, e);
            }

            LOGGER.info(" Loaded {} linx plots from {}", linxDriverPlots.size(), linxPlotDir);
        }

        LOGGER.info("Loading plots");

        String referenceBqrPlot = referenceSampleId == null
                ? null
                : processVisualisation(
                        String.format("%s/redux/%s.redux.bqr.png", pipelineOutputDir, referenceSampleId));

        String tumorBqrPlot = processVisualisation(String.format("%s/redux/%s.redux.bqr.png", pipelineOutputDir, tumorSampleId));

        String purplePlotBasePath = String.format("%s/purple/plot/%s", pipelineOutputDir, tumorSampleId);
        String purpleInputPlot = processVisualisation(purplePlotBasePath + ".input.png");
        String purpleFinalCircosPlot = processVisualisation(purplePlotBasePath + ".circos.png");
        String purpleClonalityPlot = processVisualisation(purplePlotBasePath + ".somatic.clonality.png");
        String purpleCopyNumberPlot = processVisualisation(purplePlotBasePath + ".copynumber.png");
        String purpleVariantCopyNumberPlot = processVisualisation(purplePlotBasePath + ".somatic.png");
        String purplePurityRangePlot = processVisualisation(purplePlotBasePath + ".purity.range.png");
        String purpleKataegisPlot = processVisualisation(purplePlotBasePath + ".somatic.rainfall.png");

        String cuppaSummaryPlot = hasCuppa ?
                processVisualisation(String.format("%s/cuppa/%s.cuppa.vis.png", pipelineOutputDir, tumorSampleId))
                : null;

        return VisualisationFilesBuilder.builder()
                .referenceBqrPlot(referenceBqrPlot)
                .tumorBqrPlot(tumorBqrPlot)
                .purpleInputPlot(purpleInputPlot)
                .purpleFinalCircosPlot(purpleFinalCircosPlot)
                .purpleClonalityPlot(purpleClonalityPlot)
                .purpleCopyNumberPlot(purpleCopyNumberPlot)
                .purpleVariantCopyNumberPlot(purpleVariantCopyNumberPlot)
                .purplePurityRangePlot(purplePurityRangePlot)
                .purpleKataegisPlot(purpleKataegisPlot)
                .qseePlot(processVisualisation(String.format("%s/qsee/%s.qsee.vis.report.pdf", pipelineOutputDir, tumorSampleId)))
                .linxDriverPlots(linxDriverPlots)
                .cuppaSummaryPlot(cuppaSummaryPlot)
                .build();
    }

    public String processVisualisation(final Path path)
    {
        if (!path.toFile().exists())
        {
            throw new IllegalArgumentException("Required plot file does not exist: " + path);
        }
        
        // copy the file to the plot directory
        try
        {
            LOGGER.info("Copying plot file: {}", path.getFileName());
            Path targetPath = Path.of(visDir + "/" + path.getFileName());
            Files.copy(path, targetPath, StandardCopyOption.REPLACE_EXISTING);
        }
        catch (IOException e)
        {
            throw new UncheckedIOException(e);
        }

        return path.getFileName().toString();
    }

    public String processVisualisation(final String path)
    {
        return processVisualisation(Path.of(path));
    }
}