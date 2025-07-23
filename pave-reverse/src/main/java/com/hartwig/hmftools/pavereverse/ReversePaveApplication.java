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

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class ReversePaveApplication
{
    private final ReversePaveConfig mConfig;
    private final ReversePave reversePave;

    public ReversePaveApplication(ConfigBuilder configBuilder)
    {
        mConfig = new ReversePaveConfig(configBuilder);
        reversePave = new ReversePave(mConfig.EnsemblCache, mConfig.EnsembleDataDir, mConfig.RefGenome);
    }

    public void run()
    {
        if(mConfig.Mode.equals(ROUND_TRIP_MODE))
        {
            roundTrip();
        }
        else if(mConfig.Mode.equals(SERVE_JSON_MODE))
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
        ProcessServeData processServeData = new ProcessServeData(reversePave, mConfig.RefGenVersion);
        try
        {
            processServeData.checkServeData(mConfig.ServeJsonInputFile, mConfig.TsvOuputFile);
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
            batchProcessor.process(mConfig.TsvInputFile, mConfig.TsvOuputFile);
        }
        catch(IOException e)
        {
            RPV_LOGGER.error("Failed to process batch", e);
        }
    }

    private void roundTrip()
    {
        RoundTripChecker checker = new RoundTripChecker(reversePave);

        VcfFileReader vcfFileReader = new VcfFileReader(mConfig.VcfFile, true);
        VCFHeader inputHeader = vcfFileReader.vcfHeader();

        for(VariantContext context : vcfFileReader.iterator())
        {
            VariantContext fullContext = context.fullyDecode(inputHeader, false);
            VariantImpact variantImpact = VariantImpactSerialiser.fromVariantContext(fullContext);
            checker.compareActualChangesWithCalculated(fullContext, variantImpact);
        }
        vcfFileReader.close();
        checker.printResults();
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        ReversePaveConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ReversePaveApplication reversePaveApplication = new ReversePaveApplication(configBuilder);
        reversePaveApplication.run();
    }
}

