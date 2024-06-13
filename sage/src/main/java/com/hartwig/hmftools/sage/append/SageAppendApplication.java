package com.hartwig.hmftools.sage.append;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.appendHeader;
import static com.hartwig.hmftools.sage.vcf.VcfTags.VERSION_META_DATA;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.sage.SageCommon;
import com.hartwig.hmftools.sage.evidence.FragmentLengths;
import com.hartwig.hmftools.sage.pipeline.ChromosomePartition;
import com.hartwig.hmftools.sage.bqr.BaseQualityRecalibration;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

public class SageAppendApplication
{
    private final SageAppendConfig mConfig;
    private final IndexedFastaSequenceFile mRefGenome;
    private final FragmentLengths mFragmentLengths;

    private static final double MIN_PRIOR_VERSION = 2.8;

    public SageAppendApplication(final ConfigBuilder configBuilder)
    {
        final VersionInfo version = fromAppName(APP_NAME);
        mConfig = new SageAppendConfig(version.version(), configBuilder);
        mFragmentLengths = new FragmentLengths(mConfig.Common);

        if(!mConfig.Common.isValid())
        {
            SG_LOGGER.error("invalid config");
            System.exit(1);
        }

        IndexedFastaSequenceFile refFastaSeqFile = null;

        try
        {
            refFastaSeqFile = new IndexedFastaSequenceFile(new File(mConfig.Common.RefGenomeFile));
        }
        catch (IOException e)
        {
            SG_LOGGER.error("Reference file loading failed: {}", e.toString());
            System.exit(1);
        }

        mRefGenome = refFastaSeqFile;
    }

    public void run() throws IOException, ExecutionException, InterruptedException
    {
        if(mConfig.InputVcf.equals(mConfig.Common.OutputFile))
        {
            SG_LOGGER.error("input and output VCFs must be different");
            System.exit(1);
        }

        if(mConfig.Common.ReferenceIds.isEmpty())
        {
            SG_LOGGER.error("missing reference Id must be supplied");
            System.exit(1);
        }

        SG_LOGGER.info("appending variants with {} sample BAM(s)", mConfig.Common.ReferenceBams.size());

        SG_LOGGER.info("reading and validating file: {}", mConfig.InputVcf);

        long startTimeMs = System.currentTimeMillis();

        VcfFileReader vcfFileReader = new VcfFileReader(mConfig.InputVcf);

        if(!vcfFileReader.fileValid())
        {
            SG_LOGGER.error("invalid input VCF({})", mConfig.InputVcf);
            System.exit(1);
        }

        VCFHeader inputHeader = vcfFileReader.vcfHeader();

        if(!validateInputHeader(inputHeader))
        {
            System.exit(1);
        }

        final List<VariantContext> existingVariants = Lists.newArrayList();

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            VariantContext variant = variantContext.fullyDecode(inputHeader, false);

            if(mConfig.FilterToGenes)
            {
                VariantImpact variantImpact = VariantImpactSerialiser.fromVariantContext(variant);

                if(variantImpact == null || variantImpact.GeneName.isEmpty())
                    continue;
            }

            if(!mConfig.Common.SpecificPositions.isEmpty())
            {
                if(mConfig.Common.SpecificPositions.stream().noneMatch(x -> x.matches(variant.getContig(), variant.getStart())))
                {
                    continue;
                }
            }

            existingVariants.add(variant);
        }

        vcfFileReader.close();

        SG_LOGGER.info("loaded {} variants", existingVariants.size());

        SG_LOGGER.info("writing to file: {}", mConfig.Common.OutputFile);
        VariantVCF outputVCF = new VariantVCF(mRefGenome, mConfig.Common, inputHeader);

        if(existingVariants.isEmpty())
        {
            outputVCF.close();
            SG_LOGGER.info("writing empty output VCF", existingVariants.size());
            return;
        }

        SageCommon.setReadLength(mConfig.Common, Collections.emptyMap(), mConfig.Common.ReferenceBams.get(0));

        BaseQualityRecalibration baseQualityRecalibration = new BaseQualityRecalibration(
                mConfig.Common, mRefGenome, "", Collections.emptyList(), Collections.emptyList());

        if(mConfig.Common.BQR.ExcludeKnown)
            baseQualityRecalibration.setKnownVariants(existingVariants);

        baseQualityRecalibration.produceRecalibrationMap();

        if(!baseQualityRecalibration.isValid())
            System.exit(1);

        final Map<String, BqrRecordMap> recalibrationMap = baseQualityRecalibration.getSampleRecalibrationMap();

        final ChromosomePartition chromosomePartition = new ChromosomePartition(mConfig.Common, mRefGenome);

        for(final SAMSequenceRecord samSequenceRecord : dictionary().getSequences())
        {
            final String chromosome = samSequenceRecord.getSequenceName();

            if(!mConfig.Common.processChromosome(chromosome))
                continue;

            SG_LOGGER.info("processing chromosome({})", chromosome);

            final List<VariantContext> chromosomeVariants = existingVariants.stream()
                    .filter(x -> x.getContig().equals(chromosome)).collect(Collectors.toList());

            List<ChrBaseRegion> chrBaseRegions = chromosomePartition.partition(chromosome);

            List<RegionAppendTask> regionTasks = Lists.newArrayList();

            for(int i = 0; i < chrBaseRegions.size(); ++i)
            {
                ChrBaseRegion region = chrBaseRegions.get(i);

                final List<VariantContext> regionVariants = chromosomeVariants.stream()
                        .filter(x -> positionWithin(x.getStart(), region.start(), region.end()))
                        .collect(Collectors.toList());

                if(regionVariants.isEmpty())
                    continue;

                regionTasks.add(new RegionAppendTask(i, region, regionVariants, mConfig, mRefGenome, recalibrationMap, mFragmentLengths));
            }

            final List<Callable> callableList = regionTasks.stream().collect(Collectors.toList());
            if(!TaskExecutor.executeTasks(callableList, mConfig.Common.Threads))
            {
                System.exit(1);
            }

            for(RegionAppendTask regionTask : regionTasks)
            {
                final List<VariantContext> updatedVariants = regionTask.finalVariants();
                updatedVariants.forEach(outputVCF::write);
            }
        }

        outputVCF.close();
        mFragmentLengths.close();

        mRefGenome.close();

        SG_LOGGER.info("SageAppend complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private boolean validateInputHeader(VCFHeader header)
    {
        double oldVersion = sageVersion(header);
        if(Doubles.lessThan(oldVersion, MIN_PRIOR_VERSION))
        {
            SG_LOGGER.error("Sage VCF version({}) older than required({})", oldVersion, MIN_PRIOR_VERSION);
            return false;
        }

        final Set<String> existingSamples = existingSamples(header);

        StringJoiner sj = new StringJoiner(", ");
        existingSamples.forEach(x -> sj.add(x));

        SG_LOGGER.info("existing VCF samples: {}", sj.toString());

        for(String refSample : mConfig.Common.ReferenceIds)
        {
            if(existingSamples.contains(refSample))
            {
                SG_LOGGER.error("config reference sample({}) already exits in input VCF", refSample);
                return false;
            }
        }

        appendHeader(header);

        return true;
    }

    private static double sageVersion(@NotNull final VCFHeader header)
    {
        VCFHeaderLine oldVersion = header.getMetaDataLine(VERSION_META_DATA);

        if(oldVersion == null)
            return 0;

        String[] versionComponents = oldVersion.getValue().split("\\.", -1);

        try
        {
            return Double.parseDouble(versionComponents[0]) + Double.parseDouble(versionComponents[1]);
        }
        catch(Exception e)
        {
            SG_LOGGER.error("failed to parse Sage version: {}", oldVersion.getValue());
            return 0;
        }
    }

    private static Set<String> existingSamples(final VCFHeader header)
    {
        return Sets.newHashSet(header.getGenotypeSamples());
    }

    private SAMSequenceDictionary dictionary() throws IOException
    {
        final String bam = mConfig.Common.ReferenceBams.get(0);

        SamReader tumorReader = SamReaderFactory.makeDefault()
                .validationStringency(mConfig.Common.BamStringency)
                .referenceSource(new ReferenceSource(mRefGenome)).open(new File(bam));

        SAMSequenceDictionary dictionary = tumorReader.getFileHeader().getSequenceDictionary();
        tumorReader.close();
        return dictionary;
    }


    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SageAppendConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SageAppendApplication application = new SageAppendApplication(configBuilder);

        try
        {
            application.run();
        }
        catch(Exception e)
        {
            SG_LOGGER.error("error: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
}
