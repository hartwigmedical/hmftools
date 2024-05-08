package com.hartwig.hmftools.pave;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.APP_NAME;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.pave.annotation.ReferenceData;

import org.jetbrains.annotations.NotNull;

public class PaveApplication
{
    private final PaveConfig mConfig;
    private final ReferenceData mReferenceData;

    private VcfWriter mVcfWriter;
    private final TranscriptWriter mTranscriptWriter;

    public PaveApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new PaveConfig(configBuilder);

        mReferenceData = new ReferenceData(mConfig, configBuilder);

        mVcfWriter = initialiseVcfWriter();

        mTranscriptWriter = new TranscriptWriter(mConfig);

        try
        {
            final VersionInfo version = fromAppName(APP_NAME);
            version.write(mConfig.OutputDir);
        }
        catch(IOException e)
        {
            System.exit(1);
        }
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            PV_LOGGER.error("invalid config, exiting");
            System.exit(1);
        }

        if(!mReferenceData.isValid())
        {
            PV_LOGGER.error("invalid reference data, exiting");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        List<ChromosomeTask> chromosomeTasks = Lists.newArrayList();
        List<String> initialRefChromosomes = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificRegions.isEmpty() && mConfig.SpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chrStr)))
            {
                mVcfWriter.onChromosomeComplete(chromosome);
                continue;
            }

            ChromosomeTask chromosomeTask = new ChromosomeTask(chromosome, mConfig, mReferenceData, mVcfWriter, mTranscriptWriter);

            chromosomeTasks.add(chromosomeTask);

            // initialise the reference data for the set of chromosomes which will be processed immediately to avoid data locking
            if(initialRefChromosomes.size() < max(mConfig.Threads, 1))
                initialRefChromosomes.add(chrStr);
        }

        // PV_LOGGER.debug("initialising reference data");
        mReferenceData.initialiseChromosomeData(initialRefChromosomes, mConfig.Threads);

        PV_LOGGER.info("sample({}) processing VCF file({})", mConfig.SampleId, mConfig.VcfFile);

        final List<Callable> callableList = chromosomeTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
        {
            System.exit(1);
        }

        mTranscriptWriter.close();
        mVcfWriter.close();

        PV_LOGGER.info("Pave complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private VcfWriter initialiseVcfWriter()
    {
        // append 'pave' to the input vcf file name if not specified
        String outputVcfFilename;

        if(mConfig.OutputVcfFile != null)
        {
            outputVcfFilename = mConfig.OutputVcfFile; // assumes includes path
        }
        else
        {
            String[] fileItems = mConfig.VcfFile.split("/");
            String filename = fileItems[fileItems.length - 1];
            int extensionIndex = filename.indexOf(".vcf");
            outputVcfFilename = mConfig.OutputDir + filename.substring(0, extensionIndex) + ".pave" + filename.substring(extensionIndex);

            if(!outputVcfFilename.endsWith(".gz")) // always writes zipped VCF even if input VCF isn't zipped
                outputVcfFilename += ".gz";
        }

        PV_LOGGER.info("writing VCF file({})", outputVcfFilename);

        VcfWriter vcfWriter = new VcfWriter(outputVcfFilename, mConfig.VcfFile);

        vcfWriter.writeHeader(mReferenceData, mConfig.SetReportable);

        return vcfWriter;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PaveConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PaveApplication paveApplication = new PaveApplication(configBuilder);
        paveApplication.run();
    }
}
