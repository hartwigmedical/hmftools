package com.hartwig.hmftools.pavereverse;

import static com.hartwig.hmftools.pavereverse.ReversePaveConfig.ROUND_TRIP_MODE;
import static com.hartwig.hmftools.pavereverse.ReversePaveConfig.RPV_LOGGER;
import static com.hartwig.hmftools.pavereverse.ReversePaveConfig.SERVE_JSON_MODE;
import static com.hartwig.hmftools.pavereverse.ReversePaveConstants.APP_NAME;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.pavereverse.batch.BatchProcessor;
import com.hartwig.hmftools.pavereverse.roundtrip.RoundTripChecker;
import com.hartwig.hmftools.pavereverse.serve.ProcessServeData;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class ReversePaveApplication
{
    private final ReversePaveConfig mConfig;
    private final ReversePave reversePave;

    public ReversePaveApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new ReversePaveConfig(configBuilder);
        reversePave = new ReversePave(mConfig.mEnsemblCache, mConfig.mEnsembleDataDir, mConfig.mRefGenome);
    }

    public void run()
    {
        if(mConfig.mode.equals(ROUND_TRIP_MODE))
        {
            roundTrip();
        }
        else if(mConfig.mode.equals(SERVE_JSON_MODE))
        {
            processServeJson();
        }
        else
        {
            processBatch();
        }
    }

    private void processServeJson()
    {
        ProcessServeData processServeData = new ProcessServeData(reversePave, mConfig.mRefGenVersion);
        try
        {
            processServeData.checkServeData(mConfig.mServeJsonInputFile, mConfig.mTsvOuputFile);
        }
        catch(IOException e)
        {
            RPV_LOGGER.error("Failed to process serve json", e);
        }
    }

    private void processBatch()
    {
        BatchProcessor batchProcessor = new BatchProcessor(reversePave);
        try
        {
            batchProcessor.process(mConfig.mTsvInputFile, mConfig.mTsvOuputFile);
        }
        catch(IOException e)
        {
            RPV_LOGGER.error("Failed to process batch", e);
        }
    }

    private void roundTrip()
    {
        RoundTripChecker checker = new RoundTripChecker(reversePave);

        VcfFileReader vcfFileReader = new VcfFileReader(mConfig.mVcfFile, true);
        VCFHeader inputHeader = vcfFileReader.vcfHeader();

        for(VariantContext context : vcfFileReader.iterator())
        {
            VariantContext fullContext = context.fullyDecode(inputHeader, false);
            VariantImpact variantImpact = VariantImpactSerialiser.fromVariantContext(fullContext);
            checker.compareActualChangesWithCalculated(fullContext, variantImpact);
        }
        vcfFileReader.close();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        ReversePaveConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ReversePaveApplication reversePaveApplication = new ReversePaveApplication(configBuilder);
        reversePaveApplication.run();
    }
}

