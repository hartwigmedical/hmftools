package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.CompletionException;
import java.util.concurrent.ExecutorService;
import java.util.function.BiFunction;
import java.util.function.Supplier;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.sage.config.BaseQualityRecalibrationConfig;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class QualityRecalibrationSupplier implements Supplier<Map<String, QualityRecalibrationMap>>
{
    private final ExecutorService executorService;
    private final IndexedFastaSequenceFile refGenome;
    private final SageConfig config;

    public QualityRecalibrationSupplier(final ExecutorService executorService, final IndexedFastaSequenceFile refGenome,
            final SageConfig config)
    {
        this.executorService = executorService;
        this.refGenome = refGenome;
        this.config = config;
    }

    @NotNull
    public Map<String, QualityRecalibrationMap> get()
    {
        BaseQualityRecalibrationConfig bqrConfig = config.baseQualityRecalibrationConfig();
        if(!bqrConfig.enabled())
        {
            return disableQualityRecalibration(config);
        }

        final Map<String, QualityRecalibrationMap> result = Maps.newHashMap();
        SG_LOGGER.info("Beginning quality recalibration");

        final QualityRecalibration qualityRecalibration = new QualityRecalibration(config, executorService, refGenome);
        final List<CompletableFuture<Void>> done = Lists.newArrayList();

        final BiFunction<String, String, CompletableFuture<Void>> processSample =
                (sample, sampleBam) -> qualityRecalibration.qualityRecalibrationRecords(sampleBam).thenAccept(records ->
                {
                    try
                    {

                        final String tsvFile = config.baseQualityRecalibrationFile(sample);
                        QualityRecalibrationFile.write(tsvFile, records);
                        result.put(sample, new QualityRecalibrationMap(records));
                        SG_LOGGER.info("Writing base quality recalibration file: {}", tsvFile);
                        if(bqrConfig.plot())
                        {
                            RExecutor.executeFromClasspath("r/baseQualityRecalibrationPlot.R", tsvFile);
                        }
                    } catch(Exception e)
                    {
                        throw new CompletionException(e);
                    }
                });

        for(int i = 0; i < config.reference().size(); i++)
        {
            done.add(processSample.apply(config.reference().get(i), config.referenceBam().get(i)));
        }

        for(int i = 0; i < config.tumor().size(); i++)
        {
            done.add(processSample.apply(config.tumor().get(i), config.tumorBam().get(i)));
        }

        // Wait for all tasks to be finished
        done.forEach(CompletableFuture::join);

        return result;
    }

    private Map<String, QualityRecalibrationMap> disableQualityRecalibration(SageConfig config)
    {
        final Map<String, QualityRecalibrationMap> result = Maps.newHashMap();

        for(String sample : config.reference())
        {
            result.put(sample, new QualityRecalibrationMap(Collections.emptyList()));
        }
        for(String sample : config.tumor())
        {
            result.put(sample, new QualityRecalibrationMap(Collections.emptyList()));

        }
        return result;
    }
}
