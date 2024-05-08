package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations.REPEAT_MASK_FILE;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;
import static com.hartwig.hmftools.common.variant.GenotypeIds.fromVcfHeader;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.APP_NAME;
import static com.hartwig.hmftools.esvee.caller.CallerConfig.addConfig;
import static com.hartwig.hmftools.esvee.caller.VariantFilters.logFilterTypeCounts;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.esvee.caller.annotation.PonCache;
import com.hartwig.hmftools.esvee.caller.annotation.RepeatMaskAnnotator;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class CallerApplication
{
    private final CallerConfig mConfig;
    public final FilterConstants mFilterConstants;

    private final PonCache mPonCache;
    private final HotspotCache mHotspotCache;
    private final VariantBuilder mVariantBuilder;
    private final VariantFilters mVariantFilters;
    private final RepeatMaskAnnotator mRepeatMaskAnnotator;
    private final TargetRegions mTargetRegions;

    private int mProcessedVariants;
    private final SvDataCache mSvDataCache;

    public CallerApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new CallerConfig(configBuilder);
        mFilterConstants = FilterConstants.from(configBuilder);

        SV_LOGGER.info("loading reference data");
        mPonCache = new PonCache(configBuilder);
        mHotspotCache = new HotspotCache(configBuilder);
        mTargetRegions = new TargetRegions(configBuilder);

        mVariantBuilder = new VariantBuilder(mHotspotCache, mTargetRegions);
        mVariantFilters = new VariantFilters(mFilterConstants, mConfig.hasReference(), mConfig.hasTumor());

        mProcessedVariants = 0;
        mSvDataCache = new SvDataCache();
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

        long startTimeMs = System.currentTimeMillis();

        processVcf(mConfig.VcfFile);

        SV_LOGGER.info("Esvee caller complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void processVcf(final String vcfFile)
    {
        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        VCFHeader vcfHeader = vcfFileReader.vcfHeader();

        SV_LOGGER.info("sample({}) processing VCF({})", mConfig.SampleId, vcfFile);

        GenotypeIds genotypeIds = fromVcfHeader(vcfHeader, mConfig.ReferenceId, mConfig.SampleId);

        if(genotypeIds.TumorOrdinal < 0 || (!mConfig.ReferenceId.isEmpty() && genotypeIds.ReferenceOrdinal < 0))
        {
            SV_LOGGER.error("missing sample names in VCF: {}", vcfHeader.getGenotypeSamples());
            System.exit(1);
        }

        mVariantBuilder.setGenotypeOrdinals(genotypeIds);

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

        SV_LOGGER.info("read VCF: processedBreakends({}) unmatched({}) complete({}) hardFiltered({})",
                mProcessedVariants, mVariantBuilder.incompleteSVs(), mSvDataCache.getSvList().size(), mVariantBuilder.hardFilteredCount());

        SV_LOGGER.info("writing output VCF files to {}", mConfig.OutputDir);

        final VersionInfo version = fromAppName(APP_NAME);

        VcfWriter writer = new VcfWriter(mConfig, vcfHeader, version.version(), genotypeIds, mSvDataCache);

        if(mSvDataCache.getSvList().isEmpty())
        {
            SV_LOGGER.info("writing empty VCF");
            writer.close();
            return;
        }

        SV_LOGGER.info("applying soft-filters and realignment");

        for(SvData var : mSvDataCache.getSvList())
        {
            if(mHotspotCache.isHotspotVariant(var))
                var.markHotspot();

            mVariantFilters.applyFilters(var);
        }

        mSvDataCache.buildBreakendMap();

        SV_LOGGER.info("deduplication of paired end single breakends");
        DuplicateFinder duplicateFinder = new DuplicateFinder(mSvDataCache);

        // duplicateFinder.findDuplicateSVs(alternatePaths);

        SV_LOGGER.debug("found {} SV duplications and {} SGL duplications",
                duplicateFinder.duplicateBreakends().size(), duplicateFinder.duplicateSglBreakends().size());

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

    public void processVariant(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        // SV_LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

        ++mProcessedVariants;

        if(mProcessedVariants > 0 && (mProcessedVariants % 100000) == 0)
        {
            SV_LOGGER.debug("sample({}) processed {} variants", mConfig.SampleId, mProcessedVariants);
        }

        SvData svData = mVariantBuilder.checkCreateVariant(variant, genotypeIds);

        if(svData == null)
            return;

        // optionally filter out by config
        if(mConfig.excludeVariant(svData))
            return;

        mSvDataCache.addSvData(svData);
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
