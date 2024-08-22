package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations.REPEAT_MASK_FILE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR_DESC;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;
import static com.hartwig.hmftools.common.variant.GenotypeIds.fromVcfHeader;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.caller.CallerConfig.addConfig;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.GERMLINE_AF_THRESHOLD;
import static com.hartwig.hmftools.esvee.caller.VariantFilters.logFilterTypeCounts;
import static com.hartwig.hmftools.esvee.common.CommonUtils.withLineProximity;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.common.FileCommon.formFragmentLengthDistFilename;

import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.esvee.caller.annotation.PonCache;
import com.hartwig.hmftools.esvee.caller.annotation.RepeatMaskAnnotator;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.prep.FragmentSizeDistribution;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class CallerApplication
{
    private final CallerConfig mConfig;
    public final FilterConstants mFilterConstants;

    private final PonCache mPonCache;
    private final HotspotCache mHotspotCache;
    private final VariantFilters mVariantFilters;
    private final RepeatMaskAnnotator mRepeatMaskAnnotator;

    private int mProcessedVariants;
    private final SvDataCache mSvDataCache;

    private final long mStartTimeMs;

    public CallerApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new CallerConfig(configBuilder);
        mFilterConstants = FilterConstants.from(configBuilder);

        mStartTimeMs = System.currentTimeMillis();

        SV_LOGGER.info("loading reference data");
        mPonCache = new PonCache(configBuilder);
        mHotspotCache = new HotspotCache(configBuilder);

        String inputDir = FileWriterUtils.pathFromFile(mConfig.VcfFile);
        String fragLengthFilename = formFragmentLengthDistFilename(inputDir, mConfig.SampleId);
        FragmentLengthBounds fragmentLengthBounds = FragmentSizeDistribution.loadFragmentLengthBounds(fragLengthFilename);

        SV_LOGGER.info("fragment length dist: {}", fragmentLengthBounds);

        mVariantFilters = new VariantFilters(mFilterConstants, fragmentLengthBounds);

        mProcessedVariants = 0;
        mSvDataCache = new SvDataCache(mConfig, new TargetRegions(configBuilder));
        mRepeatMaskAnnotator = new RepeatMaskAnnotator();

        if(configBuilder.hasValue(REPEAT_MASK_FILE))
        {
            if(!mRepeatMaskAnnotator.load(configBuilder.getValue(REPEAT_MASK_FILE), mConfig.RefGenVersion))
                System.exit(1);
        }
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            SV_LOGGER.error("invalid config, exiting");
            System.exit(1);
        }

        if(!mPonCache.hasValidData())
        {
            SV_LOGGER.error("invalid PON cache, exiting");
            System.exit(1);
        }

        processVcf(mConfig.VcfFile);

        SV_LOGGER.info("Esvee caller complete, mins({})", runTimeMinsStr(mStartTimeMs));
    }

    private void processVcf(final String vcfFile)
    {
        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        VCFHeader vcfHeader = vcfFileReader.vcfHeader();

        if(mConfig.ManualRefDepth > 0 && !vcfHeader.hasFormatLine(REF_DEPTH))
        {
            vcfHeader.addMetaDataLine(new VCFFormatHeaderLine(REF_DEPTH, 1, VCFHeaderLineType.Integer, REF_DEPTH_DESC));
            vcfHeader.addMetaDataLine(new VCFFormatHeaderLine(REF_DEPTH_PAIR, 1, VCFHeaderLineType.Integer, REF_DEPTH_PAIR_DESC));
        }

        SV_LOGGER.info("sample({}) processing VCF({})", mConfig.SampleId, vcfFile);

        GenotypeIds genotypeIds = fromVcfHeader(vcfHeader, mConfig.ReferenceId, mConfig.SampleId);

        if(genotypeIds.TumorOrdinal < 0 || (!mConfig.ReferenceId.isEmpty() && genotypeIds.ReferenceOrdinal < 0))
        {
            SV_LOGGER.error("missing sample names in VCF: {}", vcfHeader.getGenotypeSamples());
            System.exit(1);
        }

        mSvDataCache.setGenotypeOrdinals(genotypeIds);

        if(mConfig.GermlineOnly)
        {
            SV_LOGGER.info("germline mode ref({}: {}) tumor({}: {})",
                    genotypeIds.TumorOrdinal, genotypeIds.TumorId, genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId);
        }
        else if(mConfig.ReferenceId.isEmpty())
        {
            SV_LOGGER.info("tumor genotype info({}: {})", genotypeIds.TumorOrdinal, genotypeIds.TumorId);
        }
        else
        {
            SV_LOGGER.info("genotype info: ref({}: {}) tumor({}: {})",
                    genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId, genotypeIds.TumorOrdinal, genotypeIds.TumorId);
        }

        vcfFileReader.iterator().forEach(x -> processVariant(x, genotypeIds));

        if(mSvDataCache.incompleteSVs() == 0)
        {
            SV_LOGGER.info("loaded {} breakends, SVs({}) SGLs({})", mProcessedVariants, mSvDataCache.svCount(), mSvDataCache.sglCount());
        }
        else
        {
            SV_LOGGER.warn("loaded {} breakeds with unmatched({}) complete({}) hardFiltered({})",
                    mProcessedVariants, mSvDataCache.incompleteSVs(), mSvDataCache.getSvList().size(), mSvDataCache.hardFilteredCount());
        }

        SV_LOGGER.info("writing output VCF files to {}", mConfig.OutputDir);

        final VersionInfo version = fromAppName(APP_NAME);

        VcfWriter writer = new VcfWriter(mConfig, vcfHeader, version.version(), genotypeIds, mSvDataCache);

        if(mSvDataCache.getSvList().isEmpty())
        {
            SV_LOGGER.info("writing empty VCF");
            writer.close();
            return;
        }

        mSvDataCache.buildBreakendMap();

        markLikelyLineSites();

        SV_LOGGER.info("applying filters");

        for(Variant var : mSvDataCache.getSvList())
        {
            if(mHotspotCache.isHotspotVariant(var))
                var.markHotspot();

            markGermline(var);

            mVariantFilters.applyFilters(var);
        }

        Deduplication.deduplicateVariants(mSvDataCache.getBreakendMap());

        if(mPonCache.hasValidData())
        {
            SV_LOGGER.info("applying PON filters");
            mPonCache.annotateVariants(mSvDataCache.getSvList());
        }

        if(mRepeatMaskAnnotator.hasData())
        {
            SV_LOGGER.info("annotating with repeat-mask information");
            mRepeatMaskAnnotator.annotateVariants(mSvDataCache.getSvList());
        }

        writer.writeBreakends();
        writer.close();

        // summary logging
        if(SV_LOGGER.isDebugEnabled())
        {
            logFilterTypeCounts(mSvDataCache.getSvList());
        }
    }

    private void markGermline(final Variant var)
    {
        Breakend breakend = var.breakendStart();

        double maxGermlineAf = 0;
        double maxTumorAf = 0;

        for(Genotype genotype : breakend.Context.getGenotypes())
        {
            double af = breakend.calcAllelicFrequency(genotype);

            if(mConfig.ReferenceId.contains(genotype.getSampleName()))
                maxGermlineAf = max(maxGermlineAf, af);
            else
                maxTumorAf = max(maxTumorAf, af);
        }

        if(maxGermlineAf >= GERMLINE_AF_THRESHOLD * maxTumorAf)
            var.markGermline();
    }

    private void markLikelyLineSites()
    {
        for(List<Breakend> breakends : mSvDataCache.getBreakendMap().values())
        {
            for(int i = 0; i < breakends.size() - 1; ++i)
            {
                Breakend breakend = breakends.get(i);
                Breakend nextBreakend = breakends.get(i + 1);

                if(!breakend.isSgl() || !nextBreakend.isSgl()) // only mark SGLs, ignore breakends already linked
                    continue;

                if(!withLineProximity(breakend.Position, nextBreakend.Position, breakend.Orient, nextBreakend.Orient))
                    continue;

                // mark if either site is line
                if(breakend.Context.hasAttribute(LINE_SITE) || nextBreakend.Context.hasAttribute(LINE_SITE))
                {
                    breakend.setLineSiteBreakend(nextBreakend);
                    nextBreakend.setLineSiteBreakend(breakend);
                }
            }
        }
    }

    public void processVariant(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        // SV_LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

        ++mProcessedVariants;

        if(mProcessedVariants > 0 && (mProcessedVariants % 100000) == 0)
        {
            SV_LOGGER.debug("sample({}) processed {} variants", mConfig.SampleId, mProcessedVariants);
        }

        if(mConfig.ManualRefDepth > 0)
        {
            if(genotypeIds.hasReference())
            {
                if(!variant.getGenotype(genotypeIds.ReferenceId).hasExtendedAttribute(REF_DEPTH))
                {
                    variant.getGenotype(genotypeIds.ReferenceId).getExtendedAttributes().put(REF_DEPTH, mConfig.ManualRefDepth);
                    variant.getGenotype(genotypeIds.ReferenceId).getExtendedAttributes().put(REF_DEPTH_PAIR, mConfig.ManualRefDepth);
                }

                if(!variant.getGenotype(genotypeIds.TumorId).hasExtendedAttribute(REF_DEPTH))
                {
                    variant.getGenotype(genotypeIds.TumorId).getExtendedAttributes().put(REF_DEPTH, mConfig.ManualRefDepth);
                    variant.getGenotype(genotypeIds.TumorId).getExtendedAttributes().put(REF_DEPTH_PAIR, mConfig.ManualRefDepth);
                }
            }
        }

        mSvDataCache.processVariant(variant, genotypeIds);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CallerApplication caller = new CallerApplication(configBuilder);
        caller.run();
    }
}
