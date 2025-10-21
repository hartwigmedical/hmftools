package com.hartwig.hmftools.cobalt.e2e;

import static com.hartwig.hmftools.cobalt.metrics.MetricsConfig.BAM_FILE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.region.SpecificRegions.SPECIFIC_CHROMOSOMES;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.cobalt.metrics.FragmentMetrics;
import com.hartwig.hmftools.cobalt.metrics.WindowStatistics;
import com.hartwig.hmftools.cobalt.metrics.WindowStatisticsFile;

import org.apache.commons.io.FileUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class FragmentMetricsTest
{
    private File tempDir;
    private String sample;
    private File bamFile;
    private File outputDir;
    private List<WindowStatistics> results;

    @Before
    public void setup() throws Exception
    {
        sample = null;
        bamFile = null;
        tempDir = FileUtils.getTempDirectory();
        outputDir = new File(tempDir, "output");
        outputDir.mkdirs();
        FileUtils.cleanDirectory(outputDir);
        results = null;
    }

    @Test
    public void runTest() throws Exception
    {
        sample = "Bam21";
        bamFile = Path.of("src", "test", "resources", "bam", "chr21_reads.bam").toFile();
        runFragmentMetrics();

        assertEquals(4800, results.size());
        WindowStatistics ws0 = results.get(0);
        assertEquals(1, ws0.position());
        assertEquals(Double.NaN, ws0.median(), 0.0001);

        // The BAM has 8 pairs of reads at around 9,410,000-9,414,000
        WindowStatistics ws941 = results.get(941);
        assertEquals(9_410_001, ws941.position());
        assertEquals(8, ws941.count());
    }

    private void runFragmentMetrics() throws Exception
    {
        int argCount = 10;
        String[] args = new String[argCount];
        int index = 0;
        args[index++] = String.format("-%s", SAMPLE);
        args[index++] = String.format("%s", sample);
        args[index++] = String.format("-%s", BAM_FILE);
        args[index++] = String.format("%s", bamFile.getAbsolutePath());
        args[index++] = String.format("-%s", REF_GENOME_VERSION);
        args[index++] = String.format("%s", "37");
        args[index++] = String.format("-%s", OUTPUT_DIR);
        args[index++] = String.format("%s", outputDir.getAbsolutePath());
        args[index++] = String.format("-%s", SPECIFIC_CHROMOSOMES);
        args[index] = String.format("%d", 21);

        FragmentMetrics.main(args);

        String resultsFile = WindowStatisticsFile.fileName(outputDir.getAbsolutePath(), sample);
        results = WindowStatisticsFile.read(resultsFile);
    }
}