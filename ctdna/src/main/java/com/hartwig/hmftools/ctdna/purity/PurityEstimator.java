package com.hartwig.hmftools.ctdna.purity;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class PurityEstimator
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;

    public PurityEstimator(final CommandLine cmd)
    {
        mConfig = new PurityConfig(cmd);

        mResultsWriter = new ResultsWriter(mConfig);
    }

    public void run()
    {
        if(mConfig.PatientIds.isEmpty() || mConfig.PatientIds.size() != mConfig.SomaticVcfs.size())
        {
            CT_LOGGER.error("missing patient IDs or unmatched somatic VCFs in config");
            System.exit(1);
        }

        SomaticVariants somaticVariants = new SomaticVariants(mConfig, mResultsWriter);

        for(int i = 0; i < mConfig.PatientIds.size(); ++i)
        {
            String patientId = mConfig.PatientIds.get(i);
            String somaticVcf = mConfig.SomaticVcfs.get(i);

            somaticVariants.processPatientVcf(patientId, somaticVcf);
        }

        mResultsWriter.close();
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        PurityConfig.addCommandLineOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PurityEstimator purityEstimator = new PurityEstimator(cmd);
        purityEstimator.run();

        CT_LOGGER.info("Sample VCF analyser complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
