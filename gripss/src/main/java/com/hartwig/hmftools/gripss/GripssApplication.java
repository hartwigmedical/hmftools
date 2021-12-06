package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;
import static com.hartwig.hmftools.gripss.links.LinkRescue.findRescuedDsbLineInsertions;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.common.VcfUtils;
import com.hartwig.hmftools.gripss.filters.FilterConstants;
import com.hartwig.hmftools.gripss.filters.FilterType;
import com.hartwig.hmftools.gripss.filters.HotspotCache;
import com.hartwig.hmftools.gripss.filters.SoftFilters;
import com.hartwig.hmftools.gripss.links.AlternatePath;
import com.hartwig.hmftools.gripss.links.AlternatePathFinder;
import com.hartwig.hmftools.gripss.links.AssemblyLinks;
import com.hartwig.hmftools.gripss.links.DsbLinkFinder;
import com.hartwig.hmftools.gripss.links.LinkRescue;
import com.hartwig.hmftools.gripss.links.LinkStore;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
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

    private int mProcessedVariants;
    private final SvDataCache mSvDataCache;
    private final FilterCache mFilterCache;

    public GripssApplication(
            final GripssConfig config, final FilterConstants filterConstants, final RefGenomeInterface refGenome, final CommandLine cmd)
    {
        mConfig = config;
        mFilterConstants = filterConstants;

        mRefGenome = refGenome;

        GR_LOGGER.info("loading reference data");
        mPonCache = new PonCache(cmd);
        mHotspotCache = new HotspotCache(cmd);

        mVariantBuilder = new VariantBuilder(mFilterConstants, mHotspotCache);
        mSoftFilters = new SoftFilters(mFilterConstants);
        mRealigner = null;

        mProcessedVariants = 0;
        mSvDataCache = new SvDataCache();
        mFilterCache = new FilterCache();
    }

    public static GripssApplication fromCommandArgs(final CommandLine cmd)
    {
        GripssConfig config = new GripssConfig(cmd);
        FilterConstants filterConstants = FilterConstants.from(cmd);
        RefGenomeInterface refGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));

        return new GripssApplication(config, filterConstants, refGenome, cmd);
    }

    public void run()
    {
        processVcf(mConfig.VcfFile);
    }

    private void processVcf(final String vcfFile)
    {
        final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                vcfFile, new VCFCodec(), false);

        VCFHeader vcfHeader = (VCFHeader)reader.getHeader();

        GR_LOGGER.info("sample({}) processing VCF({})", mConfig.SampleId, vcfFile);

        GenotypeIds genotypeIds = VcfUtils.parseVcfSampleIds(vcfHeader, mConfig.ReferenceId, mConfig.SampleId);

        if(genotypeIds == null)
        {
            System.exit(1);
        }

        mRealigner = new BreakendRealigner(mRefGenome);

        if(!mConfig.tumorOnly())
        {
            GR_LOGGER.info("genetype info: ref({}: {}) tumor({}: {})",
                    genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId, genotypeIds.TumorOrdinal, genotypeIds.TumorId);
        }
        else
        {
            GR_LOGGER.info("tumor genetype info({}: {})", genotypeIds.TumorOrdinal, genotypeIds.TumorId);
        }

        try
        {
            reader.iterator().forEach(x -> processVariant(x, genotypeIds));
        }
        catch(IOException e)
        {
            GR_LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }

        GR_LOGGER.info("read VCF: breakends({}) unmatched({}) complete({}) hardFiltered({})",
                mProcessedVariants, mVariantBuilder.incompleteSVs(), mSvDataCache.getSvList().size(), mVariantBuilder.hardFilteredCount());

        if(mSvDataCache.getSvList().isEmpty())
            return;

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
        LinkStore assemblyLinkStore = AssemblyLinks.buildAssembledLinks(mSvDataCache.getSvList());
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

        Set<Breakend> rescuedBreakends = LinkRescue.findRescuedBreakends(dsbLinkStore, mFilterCache, false);
        rescuedBreakends.addAll(findRescuedDsbLineInsertions(dsbLinkStore, mFilterCache, mFilterConstants.MinQualRescueLine));
        rescuedBreakends.addAll(LinkRescue.findRescuedBreakends(assemblyLinkStore, mFilterCache, true));
        rescuedBreakends.addAll(LinkRescue.findRescuedBreakends(transitiveLinkStore, mFilterCache, true));

        GR_LOGGER.debug("rescued {} linked variants", rescuedBreakends.size());

        mFilterCache.updateFilters(rescuedBreakends, Sets.newHashSet());

        LinkStore combinedLinks = LinkStore.from(combinedTransitiveAssemblyLinks, dsbLinkStore);

        GR_LOGGER.info("writing output VCF files to {}", mConfig.OutputDir);

        final VersionInfo version = new VersionInfo("gripss.version");

        VcfWriter writer = new VcfWriter(mConfig, vcfHeader, version.version(), genotypeIds, mSvDataCache, mFilterCache);

        Map<Breakend,String> idPathMap = AlternatePathFinder.createPathMap(alternatePaths);

        writer.write(combinedLinks, idPathMap, vcfHeader);
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

    // testing methods only
    public void clearState()
    {
        mProcessedVariants = 0;
        mVariantBuilder.clearState();
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        GripssConfig.addCommandLineOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        GripssApplication gripss = GripssApplication.fromCommandArgs(cmd);
        gripss.run();

        GR_LOGGER.info("Gripss run complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
