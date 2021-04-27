package com.hartwig.hmftools.purple.plot;

import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.purple.config.ChartConfig;
import com.hartwig.hmftools.purple.config.PurpleConfig;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class Charts
{
    private final RCharts mRCharts;
    private final PurpleConfig mConfig;
    private final CircosCharts mCircosCharts;

    public Charts(final PurpleConfig config, final ExecutorService executorService, boolean isHg38) throws IOException
    {
        mRCharts = new RCharts(config, executorService);
        mCircosCharts = new CircosCharts(config, executorService, isHg38);

        mConfig = config;

        ChartConfig chartConfig = config.Charting;

        if(chartConfig.CircosBinary.isPresent())
        {
            createDirectory(chartConfig.CircosDirectory);
        }

        if(chartConfig.Enabled || chartConfig.CircosBinary.isPresent())
        {
            createDirectory(chartConfig.PlotDirectory);
        }
    }

    public void write(
            final String referenceId, final String sampleId, boolean plotSomatics,
            @NotNull final Gender gender, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<VariantContext> somaticVariants, @NotNull final List<StructuralVariant> structuralVariants,
            @NotNull final List<FittedRegion> regions, @NotNull final List<AmberBAF> bafs)
            throws InterruptedException, ExecutionException, IOException
    {
        final ChartConfig chartConfig = mConfig.Charting;

        mCircosCharts.write(referenceId, sampleId, gender, copyNumbers, somaticVariants, structuralVariants, regions, bafs);
        final List<Future<Integer>> futures = mCircosCharts.chartFutures();

        if(chartConfig.Enabled)
        {
            futures.addAll(mRCharts.chartFutures(sampleId, plotSomatics));
        }

        for(final Future<Integer> future : futures)
        {
            // This (intentionally) has side effect of alerting users to any exceptions
            int result = future.get();
            if(result != 0)
            {
                PPL_LOGGER.warn("Error generating charts.");
            }
        }
    }

    private void createDirectory(final String dir) throws IOException
    {
        final File output = new File(dir);
        if(!output.exists() && !output.mkdirs())
        {
            throw new IOException("Unable to create chart directory " + dir);
        }
    }
}
