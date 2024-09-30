package com.hartwig.hmftools.chord;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.r.RExecutor;


public class RunChord
{
    ChordConfig mConfig;

    private static final String SCRIPT_RESOURCE_PATH = "./extractSigPredictHRD.R";

    public RunChord(ChordConfig config)
    {
        mConfig = config;
    }

    public void run() throws IOException, InterruptedException
    {
        //Configurator.setRootLevel(Level.TRACE);

//        String outDir = "/Users/lnguyen/Desktop/chord_test/";
//        String sampleId = "COLO829v003T";
//        String snvIndVcf = "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/chord/src/test/resources/vcf/COLO829v003T.purple.somatic.vcf.gz";
//        String svVcf = "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/chord/src/test/resources/vcf/COLO829v003T.purple.sv.vcf.gz";
//        String refGenomeVsn = "HG37";

        String sampleId = mConfig.SampleIds.get(0);

        int result = RExecutor.executeFromClasspath(
                SCRIPT_RESOURCE_PATH,
                true,
                mConfig.OutputDir,
                sampleId,
                mConfig.purpleSomaticVcfFile(sampleId),
                mConfig.purpleSvVcfFile(sampleId),
                mConfig.RefGenVersion.toString()
        );

        if(result != 0)
        {
            throw new IOException("R execution failed. Failed to run CHORD");
        }
    }

    public static void main(String[] args) throws IOException, InterruptedException
    {
//        // TODO: Implement `extractSigPredictHRD.R` into this java class
//        throw new UnsupportedOperationException("Running CHORD from java is not yet implemented");

    }
}