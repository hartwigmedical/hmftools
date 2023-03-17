package com.hartwig.hmftools.ctdna.purity;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;

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
        SomaticVariants somaticVariants = new SomaticVariants(mConfig, mResultsWriter);

        for(String vcf : mConfig.SomaticVcfs)
        {
            if(!somaticVariants.processVcf(vcf))
            {
                System.exit(1);
            }
        }

        PurityContext purityContext = null;

        try
        {
            purityContext = PurityContextFile.read(mConfig.PurpleDir, mConfig.TumorId);
        }
        catch(Exception e)
        {
            CT_LOGGER.error("failed to load Purple purity: {}", e.toString());
            System.exit(1);
        }

        CopyNumberProfile copyNumberProfile = new CopyNumberProfile(mConfig);

        for(String ctDnaSample : mConfig.CtDnaSamples)
        {
            CnPurityResult cnPurityResult = copyNumberProfile.processSample(mConfig.TumorId, ctDnaSample);

            SomaticVariantResult somaticVariantResult = somaticVariants.processSample(ctDnaSample, purityContext);

            mResultsWriter.writeSampleSummary(ctDnaSample, cnPurityResult, somaticVariantResult);
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

        CT_LOGGER.info("Patient purity estimator complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
