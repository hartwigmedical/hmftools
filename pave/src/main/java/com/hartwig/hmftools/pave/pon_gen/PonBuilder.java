package com.hartwig.hmftools.pave.pon_gen;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.APP_NAME;
import static com.hartwig.hmftools.pave.pon_gen.PonConfig.MANUAL_ENTRIES;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Queue;
import java.util.StringJoiner;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.TaskQueue;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;
import com.hartwig.hmftools.pave.annotation.PonAnnotation;

import org.jetbrains.annotations.NotNull;

public class PonBuilder
{
    private final List<String> mSampleIds;
    private final PonConfig mConfig;

    private final PonAnnotation mExistingPon;
    private final ClinvarAnnotation mClinvarAnnotation;
    private final HotspotCache mHotspotCache;
    private final EnsemblDataCache mEnsemblDataCache;

    private final List<VariantPonData> mManualEntries;

    public PonBuilder(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);

        mConfig = new PonConfig(configBuilder);

        mExistingPon = new PonAnnotation(mConfig.ExistingPonFilename, false);

        mClinvarAnnotation = new ClinvarAnnotation(configBuilder);
        mClinvarAnnotation.loadData();

        mHotspotCache = new HotspotCache(configBuilder);

        if(configBuilder.hasValue(ENSEMBL_DATA_DIR))
        {
            mEnsemblDataCache = new EnsemblDataCache(configBuilder);
            mEnsemblDataCache.setRequiredData(true, false, false, true);
            mEnsemblDataCache.load(false);
        }
        else
        {
            mEnsemblDataCache = null;
        }

        mManualEntries = Lists.newArrayList();

        if(configBuilder.hasValue(MANUAL_ENTRIES))
        {
            String[] entries = configBuilder.getValue(MANUAL_ENTRIES).split(ITEM_DELIM, -1);

            for(String entry : entries)
            {
                String[] varValues = entry.split(":", 5);
                int sampleCount = Integer.parseInt(varValues[4]);

                VariantPonData variant = new VariantPonData(varValues[0], Integer.parseInt(varValues[1]), varValues[2], varValues[3]);
                variant.setSampleCount(sampleCount);

                mManualEntries.add(variant);
            }
        }
    }

    public void run()
    {
        if(mSampleIds.isEmpty() || mConfig.VcfPath == null)
        {
            PV_LOGGER.error("missing sample IDs file or VCF path in config");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        PV_LOGGER.info("generating PON from {} samples, minSamples({}) qualCutoff({}) mqfCutoff({})",
                mSampleIds.size(), mConfig.MinSamples, mConfig.QualCutoff, mConfig.MqfCutoff);

        List<String> sampleVcfs = Lists.newArrayListWithCapacity(mSampleIds.size());

        for(String sampleId : mSampleIds)
        {
            String vcfFilename = mConfig.VcfPath.replaceAll("\\*", sampleId);

            if(!Files.exists(Paths.get(vcfFilename)))
            {
                PV_LOGGER.warn("sample({}) missing VCF({})", sampleId, vcfFilename);
                continue;
            }

            sampleVcfs.add(vcfFilename);
        }

        RefGenomeCoordinates coordinates = mConfig.RefGenVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        Queue<ChrBaseRegion> regions = new ConcurrentLinkedQueue<>();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chrStr))
                continue;

            List<ChrBaseRegion> chrRegions = buildPartitions(chrStr, coordinates.length(chrStr), mConfig.PartitionSize);

            for(ChrBaseRegion region : chrRegions)
            {
                if(mConfig.SpecificChrRegions.hasFilters() && !mConfig.SpecificChrRegions.includeRegion(region))
                    continue;

                regions.add(region);
            }
        }

        List<PonThread> ponThreads = Lists.newArrayList();
        List<Thread> workers = new ArrayList<>();

        TaskQueue taskQueue = new TaskQueue(regions, "region", 100);

        for(int i = 0; i < min(regions.size(), mConfig.Threads); ++i)
        {
            PonThread ponThread = new PonThread(mConfig, sampleVcfs, taskQueue, mClinvarAnnotation, mHotspotCache, mEnsemblDataCache);
            ponThreads.add(ponThread);

            ponThread.start();
            workers.add(ponThread);
        }

        if(!runThreadTasks(workers))
            System.exit(1);

        List<VariantPonData> finalFilteredVariants = Lists.newArrayList();
        ponThreads.forEach(x -> finalFilteredVariants.addAll(x.filteredVariants()));

        writePon(finalFilteredVariants);

        PV_LOGGER.info("Pave PON building complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void writePon(final List<VariantPonData> variants)
    {
        // TODO: check for duplicates
        variants.addAll(mManualEntries);

        Collections.sort(variants, new VariantPonData.VariantSorter());

        try
        {
            String fileName = mConfig.OutputFilename != null ? mConfig.OutputFilename :
                    format("%s/somatic_pon_%d_samples.%s.tsv", mConfig.OutputDir, mSampleIds.size(), mConfig.RefGenVersion.identifier());

            PV_LOGGER.info("writing {} variants to PON file({})", variants.size(), fileName);

            BufferedWriter writer = createBufferedWriter(fileName);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME);
            sj.add(FLD_POSITION);
            sj.add(FLD_REF);
            sj.add(FLD_ALT);

            sj.add("SampleCount");
            sj.add("MaxSampleReadCount");
            sj.add("TotalReadCount");
            sj.add("SomaticHotspot");
            sj.add("GermlineHotspot");
            sj.add("ClinvarPathogenic");
            sj.add("InCodingRegion");
            sj.add("RepeatCount");

            writer.write(sj.toString());
            writer.newLine();

            for(VariantPonData variant : variants)
            {
                if(variant.sampleCount() < mConfig.MinSamples)
                    continue;

                if(mExistingPon.hasEntry(variant.Chromosome, variant.Position, variant.Ref, variant.Alt))
                    continue;

                sj = new StringJoiner(TSV_DELIM);
                sj.add(variant.Chromosome);
                sj.add(String.valueOf(variant.Position));
                sj.add(variant.Ref);
                sj.add(variant.Alt);
                sj.add(String.valueOf(variant.sampleCount()));
                sj.add(String.valueOf(variant.maxSampleReadCount()));
                sj.add(String.valueOf(variant.totalReadCount()));

                sj.add(String.valueOf(variant.isSomaticHotspot()));
                sj.add(String.valueOf(variant.isGermlineHotspot()));
                sj.add(String.valueOf(variant.isClinvarPathogenic()));
                sj.add(String.valueOf(variant.inCodingRegion()));
                sj.add(String.valueOf(variant.repeatCount()));

                writer.write(sj.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise output file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        PonConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PonBuilder ponBuilder = new PonBuilder(configBuilder);
        ponBuilder.run();
    }
}
