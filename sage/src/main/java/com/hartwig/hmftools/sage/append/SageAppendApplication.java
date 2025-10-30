package com.hartwig.hmftools.sage.append;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.config.VersionInfo.fromAppName;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LPS_APPEND_INFO;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LPS_APPEND_INFO_DESC;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConfig.isUltima;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.addGenotypeHeader;
import static com.hartwig.hmftools.sage.vcf.VcfTags.VERSION_META_DATA;

import java.io.File;
import java.io.IOException;
import java.lang.module.ModuleDescriptor.Version;
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
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.VersionInfo;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.sage.SageCommon;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.quality.BqrCache;
import com.hartwig.hmftools.sage.quality.BqrRecordMap;
import com.hartwig.hmftools.sage.evidence.FragmentLengthWriter;
import com.hartwig.hmftools.sage.pipeline.ChromosomePartition;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
import com.hartwig.hmftools.sage.seqtech.UltimaUtils;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class SageAppendApplication
{
    private final SageAppendConfig mConfig;
    private final IndexedFastaSequenceFile mRefGenome;
    private final FragmentLengthWriter mFragmentLengths;

    private static final Version MIN_PRIOR_VERSION = Version.parse("2.8");

    public SageAppendApplication(final ConfigBuilder configBuilder)
    {
        final VersionInfo version = fromAppName(APP_NAME);

        SageConfig.AppendMode = true;
        mConfig = new SageAppendConfig(version.version(), configBuilder);
        mFragmentLengths = new FragmentLengthWriter(mConfig.Common);

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
        catch(IOException e)
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

            if(!mConfig.Common.SpecificVariants.isEmpty())
            {
                if(mConfig.Common.SpecificVariants.stream().noneMatch(x -> x.matches(variant.getContig(), variant.getStart())))
                {
                    continue;
                }
            }
            else if(mConfig.Common.SpecificChrRegions.hasFilters())
            {
                if(!mConfig.Common.SpecificChrRegions.includePosition(variant.getContig(), variant.getStart()))
                    continue;
            }

            existingVariants.add(variant);
        }

        vcfFileReader.close();

        SG_LOGGER.info("loaded {} variants", existingVariants.size());

        SG_LOGGER.info("writing to file: {}", mConfig.Common.OutputFile);

        VariantVCF outputVCF = new VariantVCF(mRefGenome, mConfig.Common.ReferenceIds, inputHeader, mConfig.Common.OutputFile);

        if(existingVariants.isEmpty())
        {
            outputVCF.close();
            mFragmentLengths.close();
            SG_LOGGER.info("writing empty output VCF and fragment lengths TSV");
            return;
        }

        SageCommon.setReadLength(mConfig.Common, Collections.emptyMap(), mConfig.Common.ReferenceBams.get(0));

        BqrCache bqrCache = new BqrCache(mConfig.Common, Collections.emptyList());

        if(!bqrCache.isValid())
            System.exit(1);

        if(isUltima())
            UltimaUtils.setMaxRawQual(bqrCache.maxRawQual());

        final Map<String, BqrRecordMap> recalibrationMap = bqrCache.getSampleRecalibrationMap();

        MsiJitterCalcs msiJitterCalcs = MsiJitterCalcs.build(
                mConfig.Common.ReferenceIds, !mConfig.Common.SkipMsiJitter ? mConfig.Common.JitterBqrDir : null,
                mConfig.Common.Quality.HighDepthMode);

        ChromosomePartition chromosomePartition = new ChromosomePartition(mConfig.Common, mRefGenome);

        for(SAMSequenceRecord samSequenceRecord : dictionary().getSequences())
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

                regionTasks.add(new RegionAppendTask(
                        i, region, regionVariants, mConfig, mRefGenome, recalibrationMap, mFragmentLengths, msiJitterCalcs));
            }

            final List<Callable<Void>> callableList = regionTasks.stream().collect(Collectors.toList());
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
        Version version = sageVersion(header);

        if(version.compareTo(MIN_PRIOR_VERSION) < 0)
        {
            SG_LOGGER.error("Sage VCF version({}) older than required({})", version, MIN_PRIOR_VERSION);
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

        addGenotypeHeader(header); // called again in case new genotype fields are added and set in this version
        addAppendHeader(header);

        return true;
    }

    private static void addAppendHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFFormatHeaderLine(LPS_APPEND_INFO, 1, VCFHeaderLineType.String, LPS_APPEND_INFO_DESC));
    }


    private static Version sageVersion(@NotNull final VCFHeader header)
    {
        VCFHeaderLine version = header.getMetaDataLine(VERSION_META_DATA);

        if(version == null)
            return Version.parse("0");

        return Version.parse(version.getValue());
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
