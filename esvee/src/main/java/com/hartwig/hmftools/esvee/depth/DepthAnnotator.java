package com.hartwig.hmftools.esvee.depth;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR_DESC;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.common.FileCommon.DEPTH_VCF_SUFFIX;
import static com.hartwig.hmftools.esvee.common.FileCommon.ESVEE_FILE_ID;
import static com.hartwig.hmftools.esvee.common.FileCommon.formOutputFile;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.ExcludedRegions;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class DepthAnnotator
{
    private final DepthConfig mConfig;
    private final Map<String,Integer> mSampleVcfGenotypeIds;

    private final Map<String,List<VariantContext>> mChrVariantMap;
    private final List<ChrBaseRegion> mExcludedRegions;

    public DepthAnnotator(final ConfigBuilder configBuilder)
    {
        mConfig = new DepthConfig(configBuilder);
        mChrVariantMap = Maps.newHashMap();
        mSampleVcfGenotypeIds = Maps.newHashMap();
        mExcludedRegions = ExcludedRegions.getPolyGRegions(mConfig.RefGenVersion);
    }

    public void run()
    {
        if(mConfig.InputVcf == null || !Files.exists(Paths.get(mConfig.InputVcf)))
        {
            SV_LOGGER.error("missing VCF config or file");
            System.exit(1);
        }

        if(mConfig.BamFiles.size() != mConfig.Samples.size())
        {
            SV_LOGGER.error("inconsistent samples and BAM files");
            System.exit(1);
        }

        SV_LOGGER.info("running depth annotation for samples: {}", mConfig.Samples);

        long startTimeMs = System.currentTimeMillis();

        VcfFileReader reader = new VcfFileReader(mConfig.InputVcf);

        VCFHeader vcfHeader = reader.vcfHeader();

        if(!establishGenotypeIds(vcfHeader))
        {
            System.exit(1);
        }

        int vcfCount = 0;

        try
        {
            String currentChromosome = "";
            List<VariantContext> variantsList = null;

            for(VariantContext variantContext : reader.iterator())
            {
                ++vcfCount;

                if(excludeVariant(variantContext))
                    continue;

                String chromosome = variantContext.getContig();

                if(!chromosome.equals(currentChromosome))
                {
                    variantsList = Lists.newArrayList();
                    mChrVariantMap.put(chromosome, variantsList);
                    currentChromosome = chromosome;
                }

                // create new instances of each variant to set depth values against
                VariantContext newVariant = new VariantContextBuilder(variantContext)
                        .genotypes(variantContext.getGenotypes())
                        .filters(variantContext.getFilters())
                        .make();

                variantsList.add(newVariant);
            }
        }
        catch(Exception e)
        {
            SV_LOGGER.error("error reading vcf({}): {}", mConfig.InputVcf, e.toString());
            System.exit(1);
        }

        SV_LOGGER.info("loaded {} variants from vcf({})", vcfCount, mConfig.InputVcf);

        if(mChrVariantMap.isEmpty())
        {
            SV_LOGGER.warn("all variants filtered from vcf({})", vcfCount, mConfig.InputVcf);

            writeVcf(vcfHeader, Collections.emptyList());
            return;
        }

        if(mConfig.PerfLogTime > 0)
            analyseVariantDistribution();

        List<DepthTask> depthTasks = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            List<VariantContext> variantsList = mChrVariantMap.get(chrStr);

            if(variantsList == null)
                continue;

            DepthTask depthTask = new DepthTask(chrStr, mConfig, mSampleVcfGenotypeIds);
            depthTask.addVariants(variantsList);
            depthTasks.add(depthTask);
        }

        final List<Callable> callableList = depthTasks.stream().collect(Collectors.toList());
        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            System.exit(1);

        // write output VCF
        writeVcf(vcfHeader, depthTasks);

        SV_LOGGER.info("depth annotation complete, mins({})", runTimeMinsStr(startTimeMs));

        PerformanceCounter perfCounter = depthTasks.get(0).getPerfCounter();
        for(int i = 1; i < depthTasks.size(); ++i)
        {
            perfCounter.merge(depthTasks.get(i).getPerfCounter());
        }

        perfCounter.logStats();
    }

    private void writeVcf(final VCFHeader header, final List<DepthTask> depthTasks)
    {
        String outputVcf = formOutputFile(mConfig.OutputDir, mConfig.sampleId(), ESVEE_FILE_ID, DEPTH_VCF_SUFFIX, mConfig.OutputId);

        SV_LOGGER.info("writing VCF: {}", outputVcf);

        VariantContextWriter writer = new VariantContextWriterBuilder()
                .setReferenceDictionary(header.getSequenceDictionary())
                .setOutputFile(outputVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();

        if(!header.hasFormatLine(ALLELE_FRACTION))
            header.addMetaDataLine(new VCFFormatHeaderLine(ALLELE_FRACTION, 1, VCFHeaderLineType.Float, ALLELE_FRACTION_DESC));

        if(!header.hasFormatLine(REF_DEPTH))
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(REF_DEPTH, 1, VCFHeaderLineType.Integer, REF_DEPTH_DESC));
            header.addMetaDataLine(new VCFInfoHeaderLine(REF_DEPTH, 1, VCFHeaderLineType.Integer, REF_DEPTH_DESC));
        }
        else if(mConfig.VcfTagPrefix != null)
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(
                    mConfig.getVcfTag(REF_DEPTH), 1, VCFHeaderLineType.Integer, REF_DEPTH_DESC));
            header.addMetaDataLine(new VCFInfoHeaderLine(
                    mConfig.getVcfTag(REF_DEPTH), 1, VCFHeaderLineType.Integer, REF_DEPTH_DESC));
        }

        if(!header.hasFormatLine(REF_DEPTH_PAIR))
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(REF_DEPTH_PAIR, 1, VCFHeaderLineType.Integer, REF_DEPTH_PAIR_DESC));
            header.addMetaDataLine(new VCFInfoHeaderLine(REF_DEPTH_PAIR, 1, VCFHeaderLineType.Integer, REF_DEPTH_PAIR_DESC));
        }
        else if(mConfig.VcfTagPrefix != null)
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(
                    mConfig.getVcfTag(REF_DEPTH_PAIR), 1, VCFHeaderLineType.Integer, REF_DEPTH_PAIR_DESC));

            header.addMetaDataLine(new VCFInfoHeaderLine(
                    mConfig.getVcfTag(REF_DEPTH_PAIR), 1, VCFHeaderLineType.Integer, REF_DEPTH_PAIR_DESC));
        }

        writer.writeHeader(header);

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());
            DepthTask depthTask = depthTasks.stream().filter(x -> x.chromosome().equals(chrStr)).findFirst().orElse(null);

            if(depthTask == null)
                continue;

            depthTask.variants().forEach(x -> writer.add(x));
        }

        writer.close();
    }

    private boolean establishGenotypeIds(final VCFHeader header)
    {
        List<String> vcfSampleNames = header.getGenotypeSamples();

        for(int s = 0; s < mConfig.Samples.size(); ++s)
        {
            String sampleId = mConfig.Samples.get(s);
            boolean found = false;

            for(int i = 0; i < vcfSampleNames.size(); ++i)
            {
                String vcfSampleName = vcfSampleNames.get(i);

                if(sampleId.equals(vcfSampleName))
                {
                    mSampleVcfGenotypeIds.put(sampleId, i);
                    found = true;
                    break;
                }
            }

            if(!found)
            {
                SV_LOGGER.error("sample({}) not found in VCF genotype names", sampleId);
                return false;
            }
        }

        return true;
    }

    private void analyseVariantDistribution()
    {
        Map<Integer,Integer> groupFrequencies = Maps.newHashMap();

        int totalGroups = 0;
        int totalVariants = 0;
        long estimatedReads = 0;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            List<VariantContext> variants = mChrVariantMap.get(chrStr);

            if(variants == null)
                continue;

            int soloVariants = 0;
            int groupCount = 0;

            int index = 0;
            while(index < variants.size())
            {
                VariantContext variant = variants.get(index);

                int variantCount = 1;

                int posStart = variant.getStart();
                int posEnd = posStart;
                int nextIndex = index + 1;
                while(nextIndex < variants.size())
                {
                    VariantContext nextVariant = variants.get(nextIndex);

                    if(nextVariant.getStart() - posEnd > mConfig.ProximityDistance)
                        break;

                    posEnd = nextVariant.getStart();
                    ++variantCount;
                    ++nextIndex;
                }

                Integer countFrequency = groupFrequencies.get(variantCount);
                groupFrequencies.put(variantCount, countFrequency != null ? countFrequency + 1 : 1);
                ++groupCount;
                ++totalGroups;
                totalVariants += variantCount;

                if(variantCount == 1)
                    ++soloVariants;

                estimatedReads += (posEnd - posStart + 2000);

                index += variantCount;
            }

            SV_LOGGER.debug("chr({}) variants({}) group({}) soloVariants({} pct={})",
                    chrStr, variants.size(), groupCount, soloVariants, format("%.3f", soloVariants / (double)variants.size()));
        }

        int largeGroupCount = 0;
        int largeVariantsCount = 0;

        for(Map.Entry<Integer,Integer> entry : groupFrequencies.entrySet())
        {
            if(entry.getKey() >= 25)
            {
                ++largeGroupCount;
                largeVariantsCount += entry.getKey() * entry.getValue();
            }
            else
            {
                SV_LOGGER.debug("group count({}) frequency({}) total variants({})",
                        entry.getKey(), entry.getValue(), entry.getKey() * entry.getValue());
            }
        }

        SV_LOGGER.debug("large group count({}) total variants({})", largeGroupCount, largeVariantsCount);

        SV_LOGGER.debug("total variants({}) groups({}) estimated reads({})", totalVariants, totalGroups, estimatedReads);
    }

    private boolean excludeVariant(final VariantContext variant)
    {
        if(mExcludedRegions.stream().anyMatch(x -> x.containsPosition(variant.getStart())))
            return true;

        if(mConfig.SpecificRegions.isEmpty())
            return false;

        return mConfig.SpecificRegions.stream().noneMatch(x -> x.containsPosition(variant.getContig(), variant.getStart()));
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        DepthConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        DepthAnnotator bamSvSlicer = new DepthAnnotator(configBuilder);
        bamSvSlicer.run();
    }
}
