package com.hartwig.hmftools.amber.e2e;

import static com.hartwig.hmftools.amber.AmberConfig.LOCI_FILE;
import static com.hartwig.hmftools.amber.AmberConfig.USE_NEW_SEGMENTER;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.amber.AmberApplication;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.pcf.PcfSegment;

import org.apache.commons.io.FileUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class RunAmberTest
{
    private String TumorSample;
    private String ReferenceSample;
    private File TumorBamFile;
    private File ReferenceBamFile;
    private File OutputDir;
    private Multimap<Chromosome, AmberBAF> Results;


    @Before
    public void setup() throws IOException
    {
        TumorSample = null;
        ReferenceSample = null;
        TumorBamFile = null;
        ReferenceBamFile = null;
        var TempDir = new File("/Users/timlavers/work/junk/basura"); //Files.createTempDirectory("amber").toFile();
        OutputDir = new File(TempDir, "output");
        //noinspection ResultOfMethodCallIgnored
        OutputDir.mkdirs();
        FileUtils.cleanDirectory(OutputDir);
        Results = null;
    }

    @Test
    public void twoChromosomes() throws Exception
    {
        AmberScenario scenario = new AmberScenario("TwoChromosomes");
        runAmber(scenario, false);
        // Check that baf.tsv file.
        scenario.checkResults(Results);

        // Check the segmentation file.
        String segmentsFile = PCFFile.generateBAFFilename(OutputDir.getAbsolutePath(), TumorSample);
        Map<String, List<BaseRegion>> pcfData = PCFFile.loadChrBaseRegions(segmentsFile);
        assert pcfData != null;
        assertEquals(2, pcfData.size());
        List<BaseRegion> chr1Positions = pcfData.get(V38.versionedChromosome(_1));
        // The R program that does segmentation puts a couple of spurious segments
        // at the start of each chromosome.
        assertEquals(5, chr1Positions.size());
        assertEquals(49554, chr1Positions.get(0).start());
    }

    @Test
    public void newSegmenter() throws Exception
    {
        AmberScenario scenario = new AmberScenario("TwoChromosomes");
        runAmber(scenario, true);

        // Check the segmentation file.
        String segmentsFile = PCFFile.generateBAFFilename(OutputDir.getAbsolutePath(), TumorSample);
        ListMultimap<Chromosome, PcfSegment> pcfData = PCFFile.readPcfFile(segmentsFile);
        assertEquals(2, pcfData.keySet().size());
        List<PcfSegment> chr1Positions = pcfData.get(_1);
        assertEquals(3, chr1Positions.size());
        assertEquals(49554, chr1Positions.get(0).start());
    }

    private void runAmber(AmberScenario scenario, boolean useNewSegmenter) throws Exception
    {
        File sitesFile = scenario.createAmberLocationsFile(OutputDir);
        TumorBamFile = scenario.getTumorBamFile();
        TumorSample = scenario.getTumorSampleName();

        int argCount = 6;
        if(TumorBamFile != null)
        {
            argCount += 4;
        }
        if(ReferenceBamFile != null)
        {
            argCount += 4;
        }
        if(useNewSegmenter)
        {
            argCount += 1;
        }
        String[] args = new String[argCount];
        int index = 0;
        args[index++] = String.format("-%s", LOCI_FILE);
        args[index++] = String.format("%s", sitesFile.getAbsolutePath());
        args[index++] = String.format("-%s", OUTPUT_DIR);
        args[index++] = String.format("%s", OutputDir.getAbsolutePath());
        args[index++] = String.format("-%s", REF_GENOME_VERSION);
        args[index++] = String.format("%s", "38");
        if(TumorBamFile != null)
        {
            args[index++] = String.format("-%s", TUMOR);
            args[index++] = String.format("%s", TumorSample);
            args[index++] = String.format("-%s", TUMOR_BAM);
            args[index++] = String.format("%s", TumorBamFile.getAbsolutePath());
        }
        if(ReferenceBamFile != null)
        {
            args[index++] = String.format("-%s", REFERENCE);
            args[index++] = String.format("%s", ReferenceSample);
            args[index++] = String.format("-%s", REFERENCE_BAM);
            args[index++] = String.format("%s", ReferenceBamFile.getAbsolutePath());
        }
        if(useNewSegmenter)
        {
            args[index] = String.format("-%s", USE_NEW_SEGMENTER);
        }

        AmberApplication.main(args);

        File bafFile = new File(OutputDir, TumorSample + ".amber.baf.tsv.gz");
        Results = AmberBAFFile.read(bafFile.getAbsolutePath(), TumorSample != null);
//        else
//        {
//            ratioFile = new File(OutputDir, referenceSample + ".cobalt.ratio.tsv.gz");
//        }
//        assertTrue(ratioFile.exists());
//        assertTrue(ratioFile.isFile());
    }
}
