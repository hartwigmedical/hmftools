package com.hartwig.hmftools.cobalt.e2e;

import static com.hartwig.hmftools.cobalt.CobaltConfig.PCF_GAMMA;
import static com.hartwig.hmftools.cobalt.CobaltConfig.TARGET_REGION_NORM_FILE;
import static com.hartwig.hmftools.cobalt.utils.CobaltOutputsComparison.COMPARISON_VALUES_DIR;
import static com.hartwig.hmftools.cobalt.utils.CobaltOutputsComparison.ORIGINAL_VALUES_DIR;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.region.SpecificRegions.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import com.hartwig.hmftools.cobalt.CobaltApplication;
import com.hartwig.hmftools.cobalt.utils.CobaltOutputsComparison;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class CobaltBatchLocal
{
    String bamBase = "/Users/timlavers/work/scratch/datasets/pmhaem/bam/";
    String cobaltOutputBase = "/Users/timlavers/work/junk/rubbish/new/";
    String gcProfile =
            "/Users/timlavers/work/data/pipeline_resources/hmf_pipeline_resources.37_v2.0--3/dna/copy_number/GC_profile.1000bp.37.cnp";
    String panel = "/Users/timlavers/work/scratch/datasets/pmhaem/resources/cobalt_normalisation.pmh-panel-v1-1.37.tsv";

    @Test
    public void run() throws Exception
    {
//                List<String> samples = List.of("Sample_13927535");
        List<String> samples = List.of("Sample_13927535", "Sample_15307684", "Sample_15846501");
        ExecutorService executorService = Executors.newFixedThreadPool(10);
        List<Future<?>> futures = new ArrayList<>();
        samples.forEach(sample -> addExecutorTaskForSample(futures, executorService, sample));
        for(Future<?> future : futures)
        {
            try
            {
                future.get();
            }
            catch(Exception e)
            {
                executorService.shutdownNow();
                throw e;
            }
        }
        executorService.shutdown();
        executorService.awaitTermination(1, TimeUnit.HOURS);

        samples.forEach(this::runComparer);
    }

    private void addExecutorTaskForSample(List<Future<?>> futures, ExecutorService executorService, String sample)
    {
        futures.add(executorService.submit(() ->
        {
            try
            {
                runCobalt(sample);
            }
            catch(Exception e)
            {
                throw new RuntimeException(e);
            }
        }));
    }

    private void runCobalt(String sample)
    {
        String bamFile = bamBase + sample + ".bam";
        String outputDir = cobaltOutputBase + sample;
        String[] args = new String[14];
        int index = 0;
//        args[index++] = String.format("-%s", SPECIFIC_REGIONS);
//        args[index++] = String.format("%s", "1:2491001-2496001");
        args[index++] = String.format("-%s", TUMOR);
        args[index++] = String.format("%s", sample);
        args[index++] = String.format("-%s", TUMOR_BAM);
        args[index++] = String.format("%s", bamFile);
        args[index++] = String.format("-%s", GC_PROFILE);
        args[index++] = String.format("%s", gcProfile);
        args[index++] = String.format("-%s", PCF_GAMMA);
        args[index++] = String.format("%d", 50);
        args[index++] = String.format("-%s", REF_GENOME_VERSION);
        args[index++] = String.format("%s", "37");
        args[index++] = String.format("-%s", OUTPUT_DIR);
        args[index++] = String.format("%s", outputDir);
        args[index++] = String.format("-%s", TARGET_REGION_NORM_FILE);
        args[index] = String.format("%s", panel);
        //2491001

        try
        {
            CobaltApplication.main(args);
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    private void runComparer(String sample)
    {
        String originalValuesDir = "/Users/timlavers/work/junk/rubbish/baseline/" + sample;
        String comparisonValuesDir = cobaltOutputBase + sample;
        String outputDir = "/Users/timlavers/work/junk/outputs/";
        String[] args = new String[8];
        int index = 0;
        args[index++] = String.format("-%s", SAMPLE);
        args[index++] = String.format("%s", sample);
        args[index++] = String.format("-%s", ORIGINAL_VALUES_DIR);
        args[index++] = String.format("%s", originalValuesDir);
        args[index++] = String.format("-%s", COMPARISON_VALUES_DIR);
        args[index++] = String.format("%s", comparisonValuesDir);
        args[index++] = String.format("-%s", OUTPUT_DIR);
        args[index] = String.format("%s", outputDir);

        try
        {
            CobaltOutputsComparison.main(args);
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }
}
