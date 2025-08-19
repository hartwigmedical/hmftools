package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.VersionInfo.fromAppName;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConfig.isUltima;
import static com.hartwig.hmftools.sage.tinc.TincAnalyser.generateTincVcfFilename;
import static com.hartwig.hmftools.sage.tinc.TincConfig.callerTincConfig;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.VersionInfo;
import com.hartwig.hmftools.sage.evidence.FragmentLengthWriter;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.ChromosomePipeline;
import com.hartwig.hmftools.sage.quality.BqrCache;
import com.hartwig.hmftools.sage.quality.BqrRecordMap;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
import com.hartwig.hmftools.sage.seqtech.UltimaUtils;
import com.hartwig.hmftools.sage.tinc.TincAnalyser;
import com.hartwig.hmftools.sage.tinc.TincConfig;
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
    private final FragmentLengthWriter mFragmentLengths;
    private final TincConfig mTincConfig;

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

        mVcfWriter = new VcfWriter(mConfig.Common, mConfig.TumorIds, mConfig.Common.ReferenceIds, mRefData.RefGenome, mConfig.RunTinc);

        mFragmentLengths = new FragmentLengthWriter(mConfig.Common);

        SG_LOGGER.info("writing to file: {}", mConfig.Common.OutputFile);

        if(mConfig.RunTinc && !mConfig.Common.ReferenceIds.isEmpty())
        {
            mTincConfig = callerTincConfig(configBuilder, mConfig);
        }
        else
        {
            mTincConfig = null;
        }
    }

    private void run() throws IOException
    {
        long startTimeMs = System.currentTimeMillis();

        SageCommon.setReadLength(mConfig.Common, mRefData.PanelWithHotspots, mConfig.TumorBams.get(0));

        BqrCache bqrCache = new BqrCache(mConfig.Common, mConfig.TumorIds);

        if(!bqrCache.isValid())
            System.exit(1);

        if(isUltima())
            UltimaUtils.setMaxRawQual(bqrCache.maxRawQual());

        final Map<String, BqrRecordMap> recalibrationMap = bqrCache.getSampleRecalibrationMap();

        List<String> combinedSampleIds = Lists.newArrayList(mConfig.TumorIds);
        combinedSampleIds.addAll(mConfig.Common.ReferenceIds);

        MsiJitterCalcs msiJitterCalcs = MsiJitterCalcs.build(
                combinedSampleIds,
                !mConfig.Common.SkipMsiJitter ? mConfig.Common.JitterParamsDir : null,
                mConfig.Common.Quality.HighDepthMode);

        final SAMSequenceDictionary dictionary = dictionary();
        for(final SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
        {
            final String chromosome = samSequenceRecord.getSequenceName();

            if(!mConfig.Common.processChromosome(chromosome))
                continue;

            final ChromosomePipeline pipeline = new ChromosomePipeline(
                    chromosome, mConfig, mRefData, recalibrationMap, msiJitterCalcs, mPhaseSetCounter, mVcfWriter, mFragmentLengths);

            pipeline.process();
        }

        mFragmentLengths.close();
        mVcfWriter.close();

        if(mTincConfig != null && !mConfig.Common.ReferenceIds.isEmpty())
        {
            TincAnalyser tincAnalyser = new TincAnalyser(mTincConfig);
            tincAnalyser.run(mConfig.Common.Filter);

            String outputVcf = mTincConfig.RewriteVcf ? mConfig.Common.OutputFile : generateTincVcfFilename(mConfig.Common.OutputFile);
            tincAnalyser.writeVcf(mRefData.RefGenome, mConfig.Common.OutputFile, outputVcf);
        }

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