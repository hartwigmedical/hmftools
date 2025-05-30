package com.hartwig.hmftools.purple.plot;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.purple.ChartConfig;
import com.hartwig.hmftools.purple.DriverSourceData;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public class Charts
{
    private final RCharts mRCharts;
    private final PurpleConfig mConfig;
    private final CircosCharts mCircosCharts;

    public Charts(final PurpleConfig config, final ExecutorService executorService, boolean isHg38)
    {
        mRCharts = new RCharts(config, executorService);
        mConfig = config;
        mCircosCharts = config.Charting.CircosBinary != null ? new CircosCharts(config, executorService, isHg38) : null;
    }

    public void write(
            final String referenceId, final String sampleId, boolean plotSomatics,
            final Gender gender, final List<PurpleCopyNumber> copyNumbers,
            final List<VariantContextDecorator> somaticVariants, final List<StructuralVariant> structuralVariants,
            final List<ObservedRegion> regions, final List<AmberBAF> bafs, final List<DriverSourceData> driverSourceData) throws Exception
    {
        final ChartConfig chartConfig = mConfig.Charting;

        final List<Future<Integer>> chartFutures = Lists.newArrayList();

        if(mCircosCharts != null)
        {
            mCircosCharts.write(referenceId, sampleId, gender, copyNumbers, somaticVariants, structuralVariants, regions, bafs, driverSourceData);
            chartFutures.addAll(mCircosCharts.chartFutures());
        }

        if(!chartConfig.Disabled)
        {
            chartFutures.addAll(mRCharts.chartFutures(sampleId, plotSomatics));
        }

        for(final Future<Integer> future : chartFutures)
        {
            int result = future.get();
            if(result != 0)
            {
                throw new Exception("charting failed");
            }
        }
    }
}
