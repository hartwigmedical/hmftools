package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.appVersionFile;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.FragmentLengths;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.ChromosomePipeline;
import com.hartwig.hmftools.sage.bqr.BaseQualityRecalibration;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
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
    private final FragmentLengths mFragmentLengths;

    private SageApplication(final ConfigBuilder configBuilder)
    {
        final VersionInfo version = fromAppName(APP_NAME);
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

        mVcfWriter = new VcfWriter(mConfig.Common, mConfig.TumorIds, mConfig.Common.ReferenceIds, mRefData.RefGenome);

        mFragmentLengths = new FragmentLengths(mConfig.Common);

        SG_LOGGER.info("writing to file: {}", mConfig.Common.OutputFile);
    }

    private void run() throws IOException
    {
        long startTimeMs = System.currentTimeMillis();
        final Coverage coverage = new Coverage(mConfig.TumorIds, mRefData.CoveragePanel.values(), mConfig.Common);

        SageCommon.setReadLength(mConfig.Common, mRefData.PanelWithHotspots, mConfig.TumorBams.get(0));

        BaseQualityRecalibration baseQualityRecalibration = new BaseQualityRecalibration(
                mConfig.Common, mRefData.RefGenome, mConfig.PanelBed, mConfig.TumorIds, mConfig.TumorBams);
        baseQualityRecalibration.produceRecalibrationMap();

        if(!baseQualityRecalibration.isValid())
            System.exit(1);

        if(mConfig.Common.bqrRecordWritingOnly())
        {
            SG_LOGGER.info("exiting after BQR read writing");
            return;
        }

        final Map<String, BqrRecordMap> recalibrationMap = baseQualityRecalibration.getSampleRecalibrationMap();

        List<String> combinedSampleIds = Lists.newArrayList(mConfig.TumorIds);
        combinedSampleIds.addAll(mConfig.Common.ReferenceIds);
        MsiJitterCalcs msiJitterCalcs = MsiJitterCalcs.build(combinedSampleIds, mConfig.Common.JitterParamsDir);

        final SAMSequenceDictionary dictionary = dictionary();
        for(final SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
        {
            final String chromosome = samSequenceRecord.getSequenceName();

            if(!mConfig.Common.processChromosome(chromosome))
                continue;

            final ChromosomePipeline pipeline = new ChromosomePipeline(
                    chromosome, mConfig, mRefData, recalibrationMap, msiJitterCalcs, coverage, mPhaseSetCounter, mVcfWriter, mFragmentLengths);

            pipeline.process();
        }

        coverage.writeFiles(mConfig.Common.OutputFile);
        mFragmentLengths.close();

        SG_LOGGER.info("Sage complete, mins({})", runTimeMinsStr(startTimeMs));
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