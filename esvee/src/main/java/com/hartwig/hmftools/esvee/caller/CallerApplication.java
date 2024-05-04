package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations.REPEAT_MASK_FILE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INSALN;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.variant.GenotypeIds.fromVcfHeader;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.APP_NAME;
import static com.hartwig.hmftools.esvee.caller.CallerConfig.addConfig;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.esvee.caller.annotation.PonCache;
import com.hartwig.hmftools.esvee.caller.annotation.RepeatMaskAnnotation;
import com.hartwig.hmftools.esvee.caller.annotation.RepeatMaskAnnotator;
import com.hartwig.hmftools.esvee.common.FilterType;

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
    private final FilterCache mFilterCache;

    public CallerApplication(
            final CallerConfig config, final FilterConstants filterConstants, final RefGenomeInterface refGenome,
            final ConfigBuilder configBuilder)
    {
        mConfig = config;
        mFilterConstants = filterConstants;

        SV_LOGGER.info("loading reference data");
        mPonCache = new PonCache(configBuilder);
        mHotspotCache = new HotspotCache(configBuilder);
        mTargetRegions = new TargetRegions(configBuilder);

        mVariantBuilder = new VariantBuilder(mFilterConstants, mHotspotCache, mTargetRegions, mConfig.GermlineMode);
        mVariantFilters = new VariantFilters(mFilterConstants, mConfig.GermlineMode);

        mProcessedVariants = 0;
        mSvDataCache = new SvDataCache();
        mFilterCache = new FilterCache();
        mRepeatMaskAnnotator = new RepeatMaskAnnotator();

        if(configBuilder.hasValue(REPEAT_MASK_FILE))
        {
            if(!mRepeatMaskAnnotator.load(configBuilder.getValue(REPEAT_MASK_FILE), mConfig.RefGenVersion))
                System.exit(1);
        }
    }

    public static CallerApplication fromCommandArgs(final ConfigBuilder configBuilder)
    {
        CallerConfig config = new CallerConfig(configBuilder);
        FilterConstants filterConstants = FilterConstants.from(configBuilder);
        RefGenomeInterface refGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));

        return new CallerApplication(config, filterConstants, refGenome, configBuilder);
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

        if(mConfig.GermlineMode)
        {
            SV_LOGGER.info("genotype info: germline mode ref({}: {}) tumor({}: {})",
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

        final VersionInfo version = new VersionInfo("gripss.version");

        VcfWriter writer = new VcfWriter(mConfig, vcfHeader, version.version(), genotypeIds, mSvDataCache, mFilterCache);

        if(mSvDataCache.getSvList().isEmpty())
        {
            SV_LOGGER.info("writing empty VCF");
            writer.close();
            return;
        }

        SV_LOGGER.info("applying soft-filters and realignment");
        int realignedCount = 0;

        for(final SvData svData : mSvDataCache.getSvList())
        {
            mFilterCache.checkHotspotFilter(mHotspotCache, svData);

            mVariantFilters.applyFilters(svData, mFilterCache);
        }

        SV_LOGGER.info("soft-filtered({}) hotspots({}) realigned({})",
                mFilterCache.getBreakendFilters().size(), mFilterCache.getHotspots().size(), realignedCount);

        mSvDataCache.buildBreakendMap();

        SV_LOGGER.info("applying PON filters");

        for(List<Breakend> chrBreakendList : mSvDataCache.getBreakendMap().values())
        {
            for(Breakend breakend : chrBreakendList)
            {
                if(breakend == breakend.sv().breakendEnd()) // skip testing the same SV again
                    continue;

                mFilterCache.checkPonFilter(mPonCache, breakend.sv());
            }
        }

        SV_LOGGER.debug("pon filtered count({})", mFilterCache.ponFilteredCount());

        SV_LOGGER.info("finding assembly links");

        // TODO: need to load assembled links or not? for double-stranded breaks, mostlly rescue if still done
        LinkStore assemblyLinkStore = new LinkStore(); // AssemblyLinks.buildAssembledLinks(mSvDataCache.getBreakendMap());
        SV_LOGGER.debug("found {} assembly links", assemblyLinkStore.getBreakendLinksMap().size());

        /*

        SV_LOGGER.info("finding alternative paths and transitive links");

        LinkStore combinedTransitiveAssemblyLinks =LinkStore.from(assemblyLinkStore);
         */

        SV_LOGGER.info("deduplication of paired end single breakends");
        DuplicateFinder duplicateFinder = new DuplicateFinder(mSvDataCache, mFilterCache);

        // duplicateFinder.findDuplicateSVs(alternatePaths);

        mFilterCache.updateFilters(duplicateFinder.rescueBreakends(), duplicateFinder.duplicateBreakends());

        // duplicateFinder.findDuplicateSingles(combinedTransitiveAssemblyLinks);

        mFilterCache.updateFilters(Sets.newHashSet(), duplicateFinder.duplicateSglBreakends());

        SV_LOGGER.debug("found {} SV duplications and {} SGL duplications",
                duplicateFinder.duplicateBreakends().size(), duplicateFinder.duplicateSglBreakends().size());

        SV_LOGGER.info("finding double stranded break links");
        LinkStore dsbLinkStore = DsbLinkFinder.findBreaks(mSvDataCache, assemblyLinkStore, mFilterCache.getDuplicateBreakends());

        SV_LOGGER.debug("found {} double stranded breaks", dsbLinkStore.getBreakendLinksMap().size());

        SV_LOGGER.info("rescuing linked variants");

        /*
        LinkRescue linkRescue = new LinkRescue();
        linkRescue.findRescuedBreakends(dsbLinkStore, mFilterCache, false);
        linkRescue.findRescuedDsbLineInsertions(dsbLinkStore, mFilterCache, mFilterConstants.MinQualRescueLine);
        linkRescue.findRescuedBreakends(assemblyLinkStore, mFilterCache, true);
        linkRescue.findRescuedBreakends(transitiveLinkStore, mFilterCache, true);
        Set<Breakend> rescuedBreakends = linkRescue.getRescueInfo().keySet();

        SV_LOGGER.debug("rescued {} linked variants", rescuedBreakends.size());

        mFilterCache.updateFilters(rescuedBreakends, Sets.newHashSet());
        */

        if(mRepeatMaskAnnotator.hasData())
        {
            int annotated = 0;
            for(List<Breakend> chrBreakendList : mSvDataCache.getBreakendMap().values())
            {
                for(Breakend breakend : chrBreakendList)
                {
                    if(!breakend.IsStart)
                        continue;

                    if(breakend.sv().insertSequence().isEmpty())
                        continue;

                    final String alignments = breakend.Context.getAttributeAsString(INSALN, "");
                    if(alignments.isEmpty())
                        continue;

                    RepeatMaskAnnotation rmAnnotation = mRepeatMaskAnnotator.annotate(breakend.sv().insertSequence(), alignments);

                    if(rmAnnotation != null)
                    {
                        breakend.sv().setRepeatMaskAnnotation(rmAnnotation);
                        ++annotated;
                    }
                }
            }

            SV_LOGGER.debug("marked {} repeat mask annotations", annotated);
        }

        // TODO: check
        LinkStore combinedLinks = dsbLinkStore; // LinkStore.from(combinedTransitiveAssemblyLinks, dsbLinkStore);

        writer.write(combinedLinks, vcfHeader);
        writer.close();

        // summary logging
        if(SV_LOGGER.isDebugEnabled())
        {
            mFilterCache.logRescuedBreakendFilters();

            Map<FilterType,Integer> filterCounts = Maps.newHashMap();

            for(List<FilterType> filters : mFilterCache.getBreakendFilters().values())
            {
                for(FilterType type : filters)
                {
                    Integer count = filterCounts.get(type);
                    filterCounts.put(type, count != null ? count + 1 : 1);
                }
            }

            for(Map.Entry<FilterType,Integer> entry : filterCounts.entrySet())
            {
                SV_LOGGER.debug("soft filter {}: count({})", entry.getKey(), entry.getValue());
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

        CallerApplication caller = CallerApplication.fromCommandArgs(configBuilder);
        caller.run();
    }
}
