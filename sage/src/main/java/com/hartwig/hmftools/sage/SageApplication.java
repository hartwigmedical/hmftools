package com.hartwig.hmftools.sage;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageCommon.logMemoryUsage;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import com.hartwig.hmftools.common.utils.MemoryCalcs;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.ChromosomePipeline;
import com.hartwig.hmftools.sage.quality.BaseQualityRecalibration;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.vcf.VcfWriter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class SageApplication implements AutoCloseable
{
    private final SageCallConfig mConfig;
    private final ReferenceData mRefData;

    private final PhaseSetCounter mPhaseSetCounter;
    private final VcfWriter mVcfWriter;

    private SageApplication(final ConfigBuilder configBuilder)
    {
        final VersionInfo version = new VersionInfo("sage.version");
        mConfig = new SageCallConfig(version.version(), configBuilder);

        if(!mConfig.isValid())
        {
            System.exit(1);
            SG_LOGGER.error("invalid config, exiting");
        }

        mRefData = new ReferenceData(mConfig, configBuilder);

        if(!mRefData.load())
        {
            System.exit(1);
            SG_LOGGER.error("invalid reference data, exiting");
        }

        mPhaseSetCounter = new PhaseSetCounter();

        mVcfWriter = new VcfWriter(mConfig.Common.Version, mConfig.Common.OutputFile, mConfig.TumorIds, mConfig.Common.ReferenceIds, mRefData);

        SG_LOGGER.info("writing to file: {}", mConfig.Common.OutputFile);
    }

    private void run() throws IOException
    {
        long startTimeMs = System.currentTimeMillis();
        final Coverage coverage = new Coverage(mConfig.TumorIds, mRefData.CoveragePanel.values(), mConfig.Common);

        BaseQualityRecalibration baseQualityRecalibration = new BaseQualityRecalibration(
                mConfig.Common, mRefData.RefGenome, mConfig.PanelBed, mConfig.TumorIds, mConfig.TumorBams);
        baseQualityRecalibration.produceRecalibrationMap();

        if(!baseQualityRecalibration.isValid())
            System.exit(1);

        final Map<String,QualityRecalibrationMap> recalibrationMap = baseQualityRecalibration.getSampleRecalibrationMap();

        int initMemory = MemoryCalcs.calcMemoryUsage();
        logMemoryUsage(mConfig.Common.PerfWarnTime, "BQR", initMemory);
        System.gc();

        int maxTaskMemory = 0;

        final SAMSequenceDictionary dictionary = dictionary();
        for(final SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
        {
            final String chromosome = samSequenceRecord.getSequenceName();

            if(!mConfig.Common.processChromosome(chromosome))
                continue;

            final ChromosomePipeline pipeline = new ChromosomePipeline(
                    chromosome, mConfig, mRefData, recalibrationMap, coverage, mPhaseSetCounter, mVcfWriter);

            pipeline.process();
            maxTaskMemory = max(pipeline.maxMemoryUsage(), maxTaskMemory);
            System.gc();
        }

        coverage.writeFiles(mConfig.Common.OutputFile);

        SG_LOGGER.info("Sage complete, mins({})", runTimeMinsStr(startTimeMs));
        SG_LOGGER.debug("Sage memory init({}mb) max({}mb)", initMemory, maxTaskMemory);
    }

    private SAMSequenceDictionary dictionary() throws IOException
    {
        final String bam = mConfig.Common.ReferenceBams.isEmpty() ? mConfig.TumorBams.get(0) : mConfig.Common.ReferenceBams.get(0);

        SamReader tumorReader = SamReaderFactory.makeDefault()
                .validationStringency(mConfig.Common.BamStringency)
                .referenceSource(new ReferenceSource(mRefData.RefGenome))
                .open(new File(bam));

        SAMSequenceDictionary dictionary = tumorReader.getFileHeader().getSequenceDictionary();

        tumorReader.close();

        return dictionary;
    }

    @Override
    public void close() throws IOException
    {
        mRefData.RefGenome.close();
        mVcfWriter.close();
    }

    public static void main(final String... args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SageCallConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SageApplication application = new SageApplication(configBuilder);
        application.run();
        application.close();
    }
}