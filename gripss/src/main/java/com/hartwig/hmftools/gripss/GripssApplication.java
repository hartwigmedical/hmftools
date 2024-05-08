package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations.REPEAT_MASK_FILE;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.BEALN;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.gripss.GripssConfig.APP_NAME;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;
import static com.hartwig.hmftools.common.variant.GenotypeIds.fromVcfHeader;
import static com.hartwig.hmftools.gripss.GripssConfig.addConfig;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterConstants;
import com.hartwig.hmftools.gripss.filters.FilterType;
import com.hartwig.hmftools.gripss.filters.HotspotCache;
import com.hartwig.hmftools.gripss.filters.SoftFilters;
import com.hartwig.hmftools.gripss.filters.TargetRegions;
import com.hartwig.hmftools.gripss.links.AlternatePath;
import com.hartwig.hmftools.gripss.links.AlternatePathFinder;
import com.hartwig.hmftools.gripss.links.AssemblyLinks;
import com.hartwig.hmftools.gripss.links.DsbLinkFinder;
import com.hartwig.hmftools.gripss.links.LinkRescue;
import com.hartwig.hmftools.gripss.links.LinkStore;
import com.hartwig.hmftools.gripss.pon.PonCache;
import com.hartwig.hmftools.gripss.rm.RepeatMaskAnnotation;
import com.hartwig.hmftools.gripss.rm.RepeatMaskAnnotator;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class GripssApplication
{
    private final GripssConfig mConfig;
    public final FilterConstants mFilterConstants;

    private final PonCache mPonCache;
    private final HotspotCache mHotspotCache;
    private final VariantBuilder mVariantBuilder;
    private final SoftFilters mSoftFilters;
    private final RefGenomeInterface mRefGenome;
    private BreakendRealigner mRealigner;
    private final RepeatMaskAnnotator mRepeatMaskAnnotator;
    private final TargetRegions mTargetRegions;

    private int mProcessedVariants;
    private final SvDataCache mSvDataCache;
    private final FilterCache mFilterCache;

    public GripssApplication(
            final GripssConfig config, final FilterConstants filterConstants, final RefGenomeInterface refGenome,
            final ConfigBuilder configBuilder)
    {
        mConfig = config;
        mFilterConstants = filterConstants;

        mRefGenome = refGenome;

        GR_LOGGER.info("loading reference data");
        mPonCache = new PonCache(configBuilder);
        mHotspotCache = new HotspotCache(configBuilder);
        mTargetRegions = new TargetRegions(configBuilder);

        mVariantBuilder = new VariantBuilder(mFilterConstants, mHotspotCache, mTargetRegions, mConfig.GermlineMode);
        mSoftFilters = new SoftFilters(mFilterConstants, mConfig.GermlineMode);
        mRealigner = null;

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

    public static GripssApplication fromCommandArgs(final ConfigBuilder configBuilder)
    {
        GripssConfig config = new GripssConfig(configBuilder);
        FilterConstants filterConstants = FilterConstants.from(configBuilder);
        RefGenomeInterface refGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));

        return new GripssApplication(config, filterConstants, refGenome, configBuilder);
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            GR_LOGGER.error("invalid config, exiting");
            System.exit(1);
        }

        if(!mPonCache.hasValidData())
        {
            GR_LOGGER.error("invalid PON cache, exiting");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        processVcf(mConfig.VcfFile);

        GR_LOGGER.info("Gripss complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void processVcf(final String vcfFile)
    {
        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        VCFHeader vcfHeader = vcfFileReader.vcfHeader();

        GR_LOGGER.info("sample({}) processing VCF({})", mConfig.SampleId, vcfFile);

        GenotypeIds genotypeIds = fromVcfHeader(vcfHeader, mConfig.ReferenceId, mConfig.SampleId);

        if(genotypeIds.TumorOrdinal < 0 || (!mConfig.ReferenceId.isEmpty() && genotypeIds.ReferenceOrdinal < 0))
        {
            GripssConfig.GR_LOGGER.error("missing sample names in VCF: {}", vcfHeader.getGenotypeSamples());
            System.exit(1);
        }

        mVariantBuilder.setGenotypeOrdinals(genotypeIds);

        mRealigner = new BreakendRealigner(mRefGenome);

        if(mConfig.GermlineMode)
        {
            GR_LOGGER.info("genotype info: germline mode ref({}: {}) tumor({}: {})",
                    genotypeIds.TumorOrdinal, genotypeIds.TumorId, genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId);
        }
        else if(mConfig.ReferenceId.isEmpty())
        {
            GR_LOGGER.info("tumor genotype info({}: {})", genotypeIds.TumorOrdinal, genotypeIds.TumorId);
        }
        else
        {
            GR_LOGGER.info("genotype info: ref({}: {}) tumor({}: {})",
                    genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId, genotypeIds.TumorOrdinal, genotypeIds.TumorId);
        }

        vcfFileReader.iterator().forEach(x -> processVariant(x, genotypeIds));

        GR_LOGGER.info("read VCF: processedBreakends({}) unmatched({}) complete({}) hardFiltered({})",
                mProcessedVariants, mVariantBuilder.incompleteSVs(), mSvDataCache.getSvList().size(), mVariantBuilder.hardFilteredCount());

        GR_LOGGER.info("writing output VCF files to {}", mConfig.OutputDir);

        final VersionInfo version = new VersionInfo("gripss.version");

        VcfWriter writer = new VcfWriter(mConfig, vcfHeader, version.version(), genotypeIds, mSvDataCache, mFilterCache);

        if(mSvDataCache.getSvList().isEmpty())
        {
            GR_LOGGER.info("writing empty VCF");
            writer.close();
            return;
        }

        GR_LOGGER.info("applying soft-filters and realignment");
        int realignedCount = 0;

        for(final SvData svData : mSvDataCache.getSvList())
        {
            // realign breakends
            final Breakend[] breakends = svData.breakends();

            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(svData.isSgl() && se == SE_END)
                    continue;

                Breakend realignedBreakend = mRealigner.realign(breakends[se], svData.isSgl(), svData.imprecise());

                if(realignedBreakend.realigned())
                {
                    ++realignedCount;
                    breakends[se] = realignedBreakend;

                    if(!svData.isSgl())
                    {
                        int otherSe = switchIndex(se);
                        Breakend realignedRemoteBreakend = mRealigner.realignRemote(breakends[otherSe], realignedBreakend);
                        breakends[otherSe] = realignedRemoteBreakend;
                    }

                    svData.onPositionsUpdated();
                }
            }

            mFilterCache.checkHotspotFilter(mHotspotCache, svData);

            mSoftFilters.applyFilters(svData, mFilterCache);
        }

        GR_LOGGER.info("soft-filtered({}) hotspots({}) realigned({})",
                mFilterCache.getBreakendFilters().size(), mFilterCache.getHotspots().size(), realignedCount);

        mSvDataCache.buildBreakendMap();

        GR_LOGGER.info("applying PON filters");

        for(List<Breakend> chrBreakendList : mSvDataCache.getBreakendMap().values())
        {
            for(Breakend breakend : chrBreakendList)
            {
                if(breakend == breakend.sv().breakendEnd()) // skip testing the same SV again
                    continue;

                mFilterCache.checkPonFilter(mPonCache, breakend.sv());
            }
        }

        GR_LOGGER.debug("pon filtered count({})", mFilterCache.ponFilteredCount());

        GR_LOGGER.info("finding assembly links");
        LinkStore assemblyLinkStore = AssemblyLinks.buildAssembledLinks(mSvDataCache.getBreakendMap());
        GR_LOGGER.debug("found {} assembly links", assemblyLinkStore.getBreakendLinksMap().size());

        GR_LOGGER.info("finding alternative paths and transitive links");

        List<AlternatePath> alternatePaths = AlternatePathFinder.findPaths(mSvDataCache, assemblyLinkStore);
        LinkStore transitiveLinkStore = AlternatePathFinder.createLinkStore(alternatePaths);

        GR_LOGGER.debug("found {} alternate paths and {} transitive links",
                alternatePaths.size(), transitiveLinkStore.getBreakendLinksMap().size());

        LinkStore combinedTransitiveAssemblyLinks = LinkStore.from(assemblyLinkStore, transitiveLinkStore);

        GR_LOGGER.info("deduplication of paired end single breakends");
        DuplicateFinder duplicateFinder = new DuplicateFinder(mSvDataCache, mFilterCache);

        duplicateFinder.findDuplicateSVs(alternatePaths);

        mFilterCache.updateFilters(duplicateFinder.rescueBreakends(), duplicateFinder.duplicateBreakends());

        duplicateFinder.findDuplicateSingles(combinedTransitiveAssemblyLinks);

        mFilterCache.updateFilters(Sets.newHashSet(), duplicateFinder.duplicateSglBreakends());

        GR_LOGGER.debug("found {} SV duplications and {} SGL duplications",
                duplicateFinder.duplicateBreakends().size(), duplicateFinder.duplicateSglBreakends().size());

        GR_LOGGER.info("finding double stranded break links");
        LinkStore dsbLinkStore = DsbLinkFinder.findBreaks(mSvDataCache, assemblyLinkStore, mFilterCache.getDuplicateBreakends());

        GR_LOGGER.debug("found {} double stranded breaks", dsbLinkStore.getBreakendLinksMap().size());

        GR_LOGGER.info("rescuing linked variants");

        LinkRescue linkRescue = new LinkRescue();
        linkRescue.findRescuedBreakends(dsbLinkStore, mFilterCache, false);
        linkRescue.findRescuedDsbLineInsertions(dsbLinkStore, mFilterCache, mFilterConstants.MinQualRescueLine);
        linkRescue.findRescuedBreakends(assemblyLinkStore, mFilterCache, true);
        linkRescue.findRescuedBreakends(transitiveLinkStore, mFilterCache, true);
        Set<Breakend> rescuedBreakends = linkRescue.getRescueInfo().keySet();

        GR_LOGGER.debug("rescued {} linked variants", rescuedBreakends.size());

        mFilterCache.updateFilters(rescuedBreakends, Sets.newHashSet());

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

                    final String alignments = breakend.Context.getAttributeAsString(BEALN, "");
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

            GR_LOGGER.debug("marked {} repeat mask annotations", annotated);
        }

        LinkStore combinedLinks = LinkStore.from(combinedTransitiveAssemblyLinks, dsbLinkStore);

        Map<Breakend,String> idPathMap = AlternatePathFinder.createPathMap(alternatePaths);

        writer.write(combinedLinks, idPathMap, vcfHeader, linkRescue);
        writer.close();

        // summary logging
        if(GR_LOGGER.isDebugEnabled())
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
                GR_LOGGER.debug("soft filter {}: count({})", entry.getKey(), entry.getValue());
            }
        }
    }

    public void processVariant(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        // GR_LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

        ++mProcessedVariants;

        if(mProcessedVariants > 0 && (mProcessedVariants % 100000) == 0)
        {
            GR_LOGGER.debug("sample({}) processed {} variants", mConfig.SampleId, mProcessedVariants);
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

        GripssApplication gripss = GripssApplication.fromCommandArgs(configBuilder);
        gripss.run();
    }
}
