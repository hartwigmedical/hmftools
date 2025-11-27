package com.hartwig.hmftools.amber.e2e;

import static com.hartwig.hmftools.amber.AmberConfig.LOCI_FILE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import com.hartwig.hmftools.amber.AmberApplication;

import org.apache.commons.io.FileUtils;
import org.junit.Before;
import org.junit.Test;

public class RunAmberTest
{
    private File TempDir;
    private String TumorSample;
    private String ReferenceSample;
    private File TumorBamFile;
    private File ReferenceBamFile;
    private File GermlineSitesFile;
    private File OutputDir;

    @Before
    public void setup() throws IOException
    {
        TumorSample = null;
        ReferenceSample = null;
        TumorBamFile = null;
        ReferenceBamFile = null;
        TempDir = Files.createTempDirectory("amber").toFile();
        OutputDir = new File(TempDir, "output");
        //noinspection ResultOfMethodCallIgnored
        OutputDir.mkdirs();
        FileUtils.cleanDirectory(OutputDir);
    }

    @Test
    public void twoChromosomes() throws Exception
    {
        AmberScenario scenario = new AmberScenario("TwoChromosomes");
        runAmber(scenario);
    }

    private void runAmber(AmberScenario scenario) throws Exception
    {
        File sitesFile = scenario.createAmberLocationsFile(OutputDir);
        TumorBamFile = scenario.getTumorBamFile();

        int argCount = 6;
        if(TumorBamFile != null)
        {
            argCount += 4;
        }
        if(ReferenceBamFile != null)
        {
            argCount += 4;
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

        AmberApplication.main(args);

//        File ratioFile;
//        if(tumorBamFile != null)
//        {
//            ratioFile = new File(OutputDir, sample + ".cobalt.ratio.tsv.gz");
//        }
//        else
//        {
//            ratioFile = new File(OutputDir, referenceSample + ".cobalt.ratio.tsv.gz");
//        }
//        assertTrue(ratioFile.exists());
//        assertTrue(ratioFile.isFile());
    }
}
