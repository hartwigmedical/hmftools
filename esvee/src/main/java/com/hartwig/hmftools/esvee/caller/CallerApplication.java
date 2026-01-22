package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.RepeatMaskAnnotations.REPEAT_MASK_FILE;
import static com.hartwig.hmftools.common.sv.SvUtils.SV_GERMLINE_AD_THRESHOLD;
import static com.hartwig.hmftools.common.sv.SvUtils.SV_GERMLINE_AF_THRESHOLD;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR_DESC;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.VersionInfo.fromAppName;
import static com.hartwig.hmftools.common.variant.GenotypeIds.fromVcfHeader;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.caller.CallerConfig.registerConfig;
import static com.hartwig.hmftools.esvee.caller.LineChecker.adjustLineSites;
import static com.hartwig.hmftools.esvee.caller.VariantFilters.logFilterTypeCounts;
import static com.hartwig.hmftools.esvee.caller.annotation.PonCache.ARTEFACT_PON_BED_SGL_FILE;
import static com.hartwig.hmftools.esvee.caller.annotation.PonCache.ARTEFACT_PON_BED_SV_FILE;
import static com.hartwig.hmftools.esvee.caller.annotation.PonCache.GERMLINE_PON_MARGIN;
import static com.hartwig.hmftools.esvee.caller.annotation.PonCache.GERMLINE_SGL_PON_MARGIN;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.common.FileCommon.formDiscordantStatsFilename;
import static com.hartwig.hmftools.esvee.common.FileCommon.formFragmentLengthDistFilename;
import static com.hartwig.hmftools.esvee.prep.types.DiscordantStats.loadDiscordantStats;

import com.google.common.annotations.VisibleForTesting;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.sv.EsveeDiscordantStats;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.VersionInfo;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.esvee.caller.annotation.PonCache;
import com.hartwig.hmftools.esvee.caller.annotation.RepeatMaskAnnotator;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.prep.FragmentSizeDistribution;
import com.hartwig.hmftools.esvee.prep.types.DiscordantStats;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class CallerApplication
{
    private final CallerConfig mConfig;
    public final FilterConstants mFilterConstants;

    private final PonCache mPonCache;
    private final PonCache mArtefactPonCache;
    private final HotspotCache mHotspotCache;
    private final VariantFilters mVariantFilters;
    private final RepeatMaskAnnotator mRepeatMaskAnnotator;
    private final boolean mTargetedPanelMode;

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

        if(!mPonCache.hasValidData())
        {
            SV_LOGGER.error("invalid PON, exiting");
            System.exit(1);
        }

        if(configBuilder.hasValue(ARTEFACT_PON_BED_SV_FILE) || configBuilder.hasValue(ARTEFACT_PON_BED_SGL_FILE))
        {
            mArtefactPonCache = new PonCache(
                    configBuilder.getInteger(GERMLINE_PON_MARGIN),
                    configBuilder.getInteger(GERMLINE_SGL_PON_MARGIN),
                    configBuilder.getValue(ARTEFACT_PON_BED_SV_FILE),
                    configBuilder.getValue(ARTEFACT_PON_BED_SGL_FILE),
                    false);

            if(!mArtefactPonCache.hasValidData())
            {
                SV_LOGGER.error("invalid artefact PON, exiting");
                System.exit(1);
            }
        }
        else
        {
            mArtefactPonCache = null;
        }

        mHotspotCache = new HotspotCache(configBuilder);

        String fragLengthFilename = formFragmentLengthDistFilename(mConfig.PrepDir, mConfig.fileSampleId(), mConfig.OutputId);
        String discStatsFilename = formDiscordantStatsFilename(mConfig.PrepDir, mConfig.fileSampleId(), mConfig.OutputId);

        if(!Files.exists(Paths.get(fragLengthFilename)) || !Files.exists(Paths.get(discStatsFilename)))
        {
            SV_LOGGER.error("missing input files: disc-stats and frag-lengths", discStatsFilename, fragLengthFilename);
            System.exit(1);
        }

        FragmentLengthBounds fragmentLengthBounds = FragmentSizeDistribution.loadFragmentLengthBounds(fragLengthFilename);

        DiscordantStats discordantStats = loadDiscordantStats(discStatsFilename);

        SV_LOGGER.info("fragment length dist: {}", fragmentLengthBounds);

        mVariantFilters = new VariantFilters(mFilterConstants, fragmentLengthBounds, discordantStats);

        mProcessedVariants = 0;

        TargetRegions targetRegions = new TargetRegions(configBuilder.getValue(TARGET_REGIONS_BED), mConfig.RefGenVersion);
        mSvDataCache = new SvDataCache(mConfig, targetRegions);
        mTargetedPanelMode = targetRegions.hasTargetRegions();

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

        SV_LOGGER.info("sample({}) processing VCF({})", mConfig.fileSampleId(), vcfFile);

        GenotypeIds genotypeIds = fromVcfHeader(vcfHeader, mConfig.ReferenceId, mConfig.TumorId);

        if((mConfig.hasTumor() && genotypeIds.TumorOrdinal < 0) || (mConfig.hasReference() && genotypeIds.ReferenceOrdinal < 0))
        {
            SV_LOGGER.error("missing sample names in VCF: {}", vcfHeader.getGenotypeSamples());
            System.exit(1);
        }

        mSvDataCache.setGenotypeOrdinals(genotypeIds);

        if(mConfig.germlineOnly())
        {
            SV_LOGGER.info("germline mode ref({}: {}) tumor({}: {})",
                    genotypeIds.TumorOrdinal, genotypeIds.TumorId, genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId);
        }
        else if(!mConfig.hasReference())
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
            SV_LOGGER.warn("loaded {} breakends with unmatched({}) complete({}) hardFiltered({})",
                    mProcessedVariants, mSvDataCache.incompleteSVs(), mSvDataCache.getSvList().size(), mSvDataCache.hardFilteredCount());
        }

        SV_LOGGER.info("writing output VCF files to {}", mConfig.OutputDir);

        final VersionInfo version = fromAppName(APP_NAME);

        VcfWriter vcfWriter = new VcfWriter(mConfig, vcfHeader, version.version(), genotypeIds, mSvDataCache);

        if(mSvDataCache.getSvList().isEmpty())
        {
            SV_LOGGER.info("writing empty VCF");
            vcfWriter.close();
            return;
        }

        mSvDataCache.buildBreakendMap();

        LineChecker.markLineSites(mSvDataCache.getBreakendMap());

        SV_LOGGER.info("applying filters");

        for(Variant var : mSvDataCache.getSvList())
        {
            if(mHotspotCache.isHotspotVariant(var))
                var.markHotspot();

            mVariantFilters.applyFilters(var);
        }

        if(mTargetedPanelMode)
            mVariantFilters.applyAdjacentFilters(mSvDataCache.getBreakendMap());

        // set germline status and final filters based on LINE
        for(Variant var : mSvDataCache.getSvList())
        {
            markGermline(var);
        }

        for(Variant var : mSvDataCache.getSvList())
        {
            adjustLineSites(var);
        }

        Deduplication.deduplicateVariants(mSvDataCache.getBreakendMap());

        if(mPonCache.hasValidData())
        {
            SV_LOGGER.info("applying PON filters");
            mPonCache.annotateVariants(mSvDataCache.getBreakendMap());
        }

        if(mArtefactPonCache != null && mArtefactPonCache.hasValidData())
        {
            SV_LOGGER.info("applying artefacts PON filters");
            mArtefactPonCache.annotateVariants(mSvDataCache.getBreakendMap());
        }

        if(mRepeatMaskAnnotator.hasData())
        {
            SV_LOGGER.info("annotating with repeat-mask information");
            mRepeatMaskAnnotator.annotateVariants(mSvDataCache.getSvList());
        }

        vcfWriter.writeBreakends();
        vcfWriter.close();

        if(mConfig.WriteBreakendTsv)
        {
            BreakendWriter breakendWriter = new BreakendWriter(mConfig);
            breakendWriter.writeBreakends(mSvDataCache);
        }

        // summary logging
        if(SV_LOGGER.isDebugEnabled())
        {
            logFilterTypeCounts(mSvDataCache.getSvList());
        }
    }

    private void markGermline(final Variant var)
    {
        if(mConfig.germlineOnly())
        {
            var.markGermline();
            return;
        }

        if(isGermline(var, mConfig.hasReference() ? mConfig.ReferenceId : null))
            var.markGermline();
    }

    @VisibleForTesting
    public static boolean isGermline(final Variant var, @Nullable final String referenceId)
    {
        Breakend breakend = var.breakendStart();

        double maxGermlineAf = 0;
        double maxTumorAf = 0;
        int germlineAd = 0;
        int tumorAd = 0;

        for(Genotype genotype : breakend.Context.getGenotypes())
        {
            double af = breakend.calcAllelicFrequency(genotype);

            if(referenceId != null && referenceId.equals(genotype.getSampleName()))
            {
                maxGermlineAf = max(maxGermlineAf, af);
                germlineAd = breakend.fragmentCount(genotype);
            }
            else
            {
                maxTumorAf = max(maxTumorAf, af);
                tumorAd = breakend.fragmentCount(genotype);
            }
        }

        if(maxGermlineAf >= SV_GERMLINE_AF_THRESHOLD * maxTumorAf)
        {
            // also check the relative fragment counts
            double adRatio = tumorAd > 0 ? germlineAd / (double)tumorAd : 1;

            if(adRatio >= SV_GERMLINE_AD_THRESHOLD)
                return true;
        }

        return false;
    }

    public void processVariant(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        // SV_LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

        ++mProcessedVariants;

        if(mProcessedVariants > 0 && (mProcessedVariants % 100000) == 0)
        {
            SV_LOGGER.debug("sample({}) processed {} variants", mConfig.TumorId, mProcessedVariants);
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
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CallerApplication caller = new CallerApplication(configBuilder);
        caller.run();
    }
}
