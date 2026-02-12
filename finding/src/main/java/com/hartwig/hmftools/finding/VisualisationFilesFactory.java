package com.hartwig.hmftools.finding;

import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.datamodel.finding.VisualisationFiles;
import com.hartwig.hmftools.datamodel.finding.VisualisationFilesBuilder;
import com.hartwig.hmftools.datamodel.orange.OrangePlots;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jspecify.annotations.Nullable;

public class VisualisationFilesFactory
{
    private static final Logger LOGGER = LogManager.getLogger(VisualisationFilesFactory.class);

    private VisualisationFilesFactory()
    {
    }

    public static VisualisationFiles create(OrangePlots orangePlots)
    {
        LOGGER.info("Loading plots");

        List<String> linxDriverPlots = orangePlots.linxDriverPlots().stream()
                .map(VisualisationFilesFactory::processFilePath)
                .toList();

        String referenceBqrPlot = processNullableFilePath(orangePlots.sageReferenceBQRPlot());

        String tumorBqrPlot = processFilePath(orangePlots.sageTumorBQRPlot());

        String purpleInputPlot = processFilePath(orangePlots.purpleInputPlot());
        String purpleFinalCircosPlot = processFilePath(orangePlots.purpleFinalCircosPlot());
        String purpleClonalityPlot = processFilePath(orangePlots.purpleClonalityPlot());
        String purpleCopyNumberPlot = processFilePath(orangePlots.purpleCopyNumberPlot());
        String purpleVariantCopyNumberPlot = processFilePath(orangePlots.purpleVariantCopyNumberPlot());
        String purplePurityRangePlot = processFilePath(orangePlots.purplePurityRangePlot());
        String purpleKataegisPlot = processFilePath(orangePlots.purpleKataegisPlot());

        String cuppaSummaryPlot = processNullableFilePath(orangePlots.cuppaSummaryPlot());

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
                .qseePlot(null)
                .linxDriverPlots(linxDriverPlots)
                .cuppaSummaryPlot(cuppaSummaryPlot)
                .build();
    }

    @Nullable
    private static String processNullableFilePath(@Nullable final String path)
    {
        if(path != null)
        {
            return processFilePath(path);
        }
        return null;
    }

    private static String processFilePath(final String path)
    {
        return Path.of(path).getFileName().toString();
    }
}