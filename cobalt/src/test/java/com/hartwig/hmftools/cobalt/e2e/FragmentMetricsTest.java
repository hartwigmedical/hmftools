package com.hartwig.hmftools.cobalt.e2e;

import static com.hartwig.hmftools.cobalt.metrics.MetricsConfig.BAM_FILE;
import static com.hartwig.hmftools.cobalt.metrics.MetricsConfig.WINDOW_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.region.SpecificRegions.SPECIFIC_CHROMOSOMES;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;

import static org.junit.Assert.assertEquals;

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
    private String sample;
    private File bamFile;
    private File outputDir;
    private List<WindowStatistics> results;

    @Before
    public void setup() throws Exception
    {
        File tempDir = FileUtils.getTempDirectory();
        outputDir = new File(tempDir, "output");
        //noinspection ResultOfMethodCallIgnored
        outputDir.mkdirs();
        FileUtils.cleanDirectory(outputDir);
        results = null;
        sample = "Bam21";
        bamFile = Path.of("src", "test", "resources", "bam", "chr21_reads.bam").toFile();
        // The BAM has these reads:
        /*
LH00144:313:22YY5KLT3:4:1229:11151:8398:CGATATG_GTGTCAT	113	21	9411208	40	143M	20	25901176	0	GAAATTGTAGTT - low quality
LH00144:313:22YY5KLT3:8:2266:50581:27661:GCCATTA_CAGTACT	99	21	9412168	60	143M	=	9412447	422	ATATCTGTGATT - count
LH00144:313:22YY5KLT3:8:2266:50581:27661:GCCATTA_CAGTACT	147	21	9412447	60	143M	=	9412168	-422	ATTATAAA - 2nd in pair
LH00144:313:22YY5KLT3:3:1201:37906:7501:ACTAGGT_GCTAACT	163	21	9412474	60	143M	=	9412871	540	AAAATAAGATACAAAA - 2nd in pair
LH00144:313:22YY5KLT3:3:1201:37906:7501:ACTAGGT_GCTAACT	83	21	9412871	60	143M	=	9412474	-540	ACATATATTGCA - count
LH00144:313:22YY5KLT3:5:1163:37592:2473:CNS_CGAATCT_TCGTGTT	99	21	9412920	60	143M	=	9413052	275	TTTAGGTCACTT - count
LH00144:313:22YY5KLT3:5:1163:37592:2473:CNS_CGAATCT_TCGTGTT	147	21	9413052	48	143M	=	9412920	-275	CTTTTAAT - 2nd in pair
LH00144:313:22YY5KLT3:3:1144:51505:12641:CAGTACT_TGTCGTA	145	21	9413464	32	143M	2	105528597	0	CCACCCTA - low quality
LH00144:313:22YY5KLT3:4:2192:28328:29407:CNS_CACTGTA_TGTCGTC	147	21	9413728	60	54S89M	=	9413729	-88	TACACTGG - 2nd in pair
LH00144:313:22YY5KLT3:4:2192:28328:29407:CNS_CACTGTA_TGTCGTC	99	21	9413729	60	89M54S	=	9413728	88	TCAAACCC - soft clipped
LH00144:313:22YY5KLT3:8:1210:47872:12241:GTCACTC_TGTCGTC	83	21	9413792	59	38S105M	=	9413793	-104	GACTGGAG - soft clipped
LH00144:313:22YY5KLT3:8:1210:47872:12241:GTCACTC_TGTCGTC	163	21	9413793	58	106M37S	=	9413792	104	TGGACTTCTGGC - soft clipped
LH00144:313:22YY5KLT3:6:2102:35105:14595:CNS_GCTAACT_CGAATCT	163	21	9413918	60	143M	=	9413918	143	CCTTGGTG - 2nd in pair
LH00144:313:22YY5KLT3:6:2102:35105:14595:CNS_GCTAACT_CGAATCT	83	21	9413918	60	143M	=	9413918	-143	CCTT - count
LH00144:313:22YY5KLT3:4:1115:9394:23674:CACTGTG_GCTGTTG	99	21	9414019	35	143M	=	9414103	227	AGTGTTTTCTCTACAC - low quality
         */
    }

    @Test
    public void defaultWindowSize() throws Exception
    {
        runFragmentMetrics(null);

        assertEquals(48000, results.size());
        WindowStatistics ws0 = results.get(0);
        assertEquals(1, ws0.position());
        assertEquals(Double.NaN, ws0.median(), 0.0001);

        WindowStatistics ws9412 = results.get(9412);
        assertEquals(9_412_001, ws9412.position());
        assertEquals(3, ws9412.count());
        assertEquals(275, ws9412.min(), 0.0001);
        assertEquals(540, ws9412.max(), 0.0001);
        assertEquals((275.0 + 422.0 + 540.0)/3.0, ws9412.mean(), 0.0001);
        assertEquals(422.0, ws9412.median(), 0.0001);

        WindowStatistics ws9413 = results.get(9413);
        assertEquals(9_413_001, ws9413.position());
        assertEquals(1, ws9413.count());
        assertEquals(143, ws9413.min(), 0.0001);
        assertEquals(143, ws9413.max(), 0.0001);
        assertEquals(143.0, ws9413.mean(), 0.0001);
        assertEquals(143.0, ws9413.median(), 0.0001);
    }

    @Test
    public void windowSize10K() throws Exception
    {
        runFragmentMetrics(10_000);

        assertEquals(4800, results.size());
        WindowStatistics ws0 = results.get(0);
        assertEquals(1, ws0.position());
        assertEquals(Double.NaN, ws0.median(), 0.0001);

        WindowStatistics ws941 = results.get(941);
        assertEquals(9_410_001, ws941.position());
        assertEquals(4, ws941.count());
        assertEquals(143, ws941.min(), 0.0001);
        assertEquals(540, ws941.max(), 0.0001);
        assertEquals((143.0 + 275.0 + 422.0 + 540.0)/4.0, ws941.mean(), 0.0001);
        assertEquals((275.0 + 422.0)/2.0, ws941.median(), 0.0001);
    }

    private void runFragmentMetrics(Integer windowSize) throws Exception
    {
        int argCount = windowSize == null ? 10 : 12;
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
        if (windowSize != null)
        {
            args[index++] = String.format("-%s", WINDOW_SIZE);
            args[index++] = String.format("%d", windowSize);
        }
        args[index++] = String.format("-%s", SPECIFIC_CHROMOSOMES);
        args[index] = String.format("%d", 21);

        FragmentMetrics.main(args);

        String resultsFile = WindowStatisticsFile.fileName(outputDir.getAbsolutePath(), sample);
        results = WindowStatisticsFile.read(resultsFile);
    }
}