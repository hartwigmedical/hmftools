package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.gripss.GermlineUtils.GM_LOGGER;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterConstants;
import com.hartwig.hmftools.gripss.filters.FilterType;
import com.hartwig.hmftools.gripss.filters.HardFilters;
import com.hartwig.hmftools.gripss.filters.SoftFilters;
import com.hartwig.hmftools.gripss.links.AlternatePath;
import com.hartwig.hmftools.gripss.links.AlternatePathFinder;
import com.hartwig.hmftools.gripss.links.AssemblyLinks;
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
    private final HardFilters mHardFilters;
    private final SoftFilters mSoftFilters;
    private final RefGenomeInterface mRefGenome;
    private BreakendRealigner mRealigner;

    private int mProcessedVariants;
    private StructuralVariantFactory mSvFactory;
    private final Set<String> mHardFilteredVcfIds;
    private final SvDataCache mSvDataCache;
    private final Map<FilterType,Integer> mFilterCounts;

    public GripssApplication(final CommandLine cmd)
    {
        mConfig = new GripssConfig(cmd);
        mFilterConstants = FilterConstants.from(cmd);

        mSvFactory = new StructuralVariantFactory(new AlwaysPassFilter());
        mHardFilteredVcfIds = Sets.newHashSet();

        mRefGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));

        GM_LOGGER.info("loading reference data");
        mPonCache = new PonCache(cmd, mConfig.RestrictedChromosomes);
        mHotspotCache = new HotspotCache(cmd);

        mHardFilters = new HardFilters(mFilterConstants, mHotspotCache);
        mSoftFilters = new SoftFilters(mFilterConstants);
        mRealigner = null;

        mProcessedVariants = 0;
        mSvDataCache = new SvDataCache();
        mFilterCounts = Maps.newHashMap();
    }

    public void run()
    {
        processVcf(mConfig.VcfFile);
    }

    private void processVcf(final String vcfFile)
    {
        try
        {
            GM_LOGGER.info("processing germline VCF({})", vcfFile);

            // TODO: only in place in run in a cohort mode, but then would be threaded anyway?
            clearSampleData();

            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    vcfFile, new VCFCodec(), false);

            GenotypeIds genotypeIds = VcfUtils.parseVcfSampleIds((VCFHeader)reader.getHeader(), mConfig.ReferenceId, mConfig.SampleId);

            if(genotypeIds == null)
            {
                System.exit(1);
            }

            mRealigner = new BreakendRealigner(mRefGenome, genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);

            if(mConfig.tumorOnly())
            {
                GM_LOGGER.info("sample({}) genetype info: ref({}: {}) tumor({}: {})",
                        genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId, genotypeIds.TumorOrdinal, genotypeIds.TumorId);
            }
            else
            {
                GM_LOGGER.info("sample({}) tumor genetype info({}: {})", genotypeIds.TumorOrdinal, genotypeIds.TumorId);
            }

            reader.iterator().forEach(x -> processVariant(x, genotypeIds));

            GM_LOGGER.info("sample({}) read VCF: unmatched({}) results({})",
                    mConfig.SampleId, mSvFactory.unmatched().size(), mSvFactory.results().size());

            if(mSvDataCache.getSvList().isEmpty())
                return;

            for(final SvData svData : mSvDataCache.getSvList())
            {
                // realign breakends
                final Breakend[] breakends = svData.breakends();

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    Breakend realignedBreakend = mRealigner.realign(breakends[se], svData.isSgl(), svData.imprecise());

                    if(realignedBreakend.realigned())
                    {
                        breakends[se] = realignedBreakend;

                        int otherSe = switchIndex(se);
                        Breakend realignedRemoteBreakend = mRealigner.realignRemote(breakends[otherSe], realignedBreakend);
                        breakends[otherSe] = realignedRemoteBreakend;
                    }
                }

                // CHECK: annotate with the PON - or can this be done at the very end only for variants written to file?
                svData.setPonCount(mPonCache.getPonCount(svData));

                mSoftFilters.applyFilters(svData);

                if(GM_LOGGER.isDebugEnabled())
                {
                    svData.getFilters().forEach(x -> registerFilter(x));

                    svData.breakendStart().getFilters().forEach(x -> registerFilter(x));

                    if(!svData.isSgl())
                        svData.breakendEnd().getFilters().forEach(x -> registerFilter(x));
                }

                // TODO maintain or ditch?
                // softFilters.forEach(x -> registerFilter(x));
            }
        }
        catch(IOException e)
        {
            GM_LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }

        mSvDataCache.buildBreakendMap();

        GM_LOGGER.info("finding assembly links");
        LinkStore assemblyLinkStore = AssemblyLinks.buildAssembledLinks(mSvDataCache.getSvList());

        GM_LOGGER.info("finding alternative paths and transitive links");
        /*
        logger.info("Finding transitive links")
        val alternatePaths: Collection<AlternatePath> = AlternatePath(assemblyLinks, variantStore)
        val alternatePathsStringsByVcfId = alternatePaths.associate { x -> Pair(x.vcfId, x.pathString()) }
        val transitiveLinks = LinkStore(alternatePaths.flatMap { x -> x.transitiveLinks() })
        val combinedTransitiveAssemblyLinks = LinkStore(assemblyLinks, transitiveLinks)
        */
        List<AlternatePath> alternatePaths = AlternatePathFinder.findPaths(mSvDataCache, assemblyLinkStore);
        Map<String,String> idPathMap = AlternatePathFinder.createIdToPathMap(alternatePaths);
        LinkStore transitiveLinkStore = AlternatePath.createLinkStore(alternatePaths);
        LinkStore combinedTransitiveAssemblyLinks = LinkStore.from(assemblyLinkStore, transitiveLinkStore);

        GM_LOGGER.info("paired break end de-duplication");

        // val dedupPair = DedupPair(initialFilters, alternatePaths, variantStore)
        // val softFiltersAfterPairedDedup = initialFilters.update(dedupPair.duplicates, dedupPair.rescue)

        GM_LOGGER.info("single break end de-duplication");

        // val dedupSingle = DedupSingle(variantStore, softFiltersAfterPairedDedup, combinedTransitiveAssemblyLinks)
        // val softFiltersAfterSingleDedup = softFiltersAfterPairedDedup.update(dedupSingle.duplicates, setOf())

        GM_LOGGER.info("Finding double stranded break links");
        // val dsbLinks = DsbLink(variantStore, assemblyLinks, softFiltersAfterSingleDedup.duplicates())

        GM_LOGGER.info("rescuing linked variants");
        /*
        val dsbRescues = LinkRescue.rescueDsb(dsbLinks, softFiltersAfterSingleDedup, variantStore).rescues
        val dsbRescueMobileElements = LinkRescue.rescueDsbMobileElementInsertion(config.filterConfig, dsbLinks, softFiltersAfterSingleDedup, variantStore).rescues
        val assemblyRescues = LinkRescue.rescueAssembly(assemblyLinks, softFiltersAfterSingleDedup, variantStore).rescues
        val transitiveRescues = LinkRescue.rescueTransitive(transitiveLinks, softFiltersAfterSingleDedup, variantStore).rescues
        val allRescues = dsbRescues + dsbRescueMobileElements + assemblyRescues + transitiveRescues
        */

        GM_LOGGER.info("writing output VCF file: {}", mConfig.OutputVcfFile);

        /*
        logger.info("Writing file: ${config.outputVcf}")
        val combinedLinks = LinkStore(combinedTransitiveAssemblyLinks, dsbLinks)
        val finalFilters: SoftFilterStore = softFiltersAfterSingleDedup.update(setOf(), allRescues)
        fileWriter.writeHeader(version.version(), fileReader.fileHeader, outputSampleNames)
        for (variant in variantStore.selectAll()) {

            val localLinkedBy = combinedLinks[variant.vcfId]
            val remoteLinkedBy = combinedLinks[variant.mateId]
            val altPath = alternatePathsStringsByVcfId[variant.vcfId]

            val filters = finalFilters.filters(variant.vcfId, variant.mateId)
            fileWriter.writeVariant(variant.context(localLinkedBy, remoteLinkedBy, altPath, hotspots.contains(variant.vcfId), filters))
        }
        */

        // summary logging
        int hardFiltered = mFilterCounts.values().stream().mapToInt(x -> x.intValue()).sum();

        GM_LOGGER.info("sample({}) read {} variants: hardFiltered({}) cached({})",
                mConfig.SampleId, mProcessedVariants, hardFiltered, mSvDataCache.getSvList().size());

        if(GM_LOGGER.isDebugEnabled())
        {
            for(Map.Entry<FilterType,Integer> entry : mFilterCounts.entrySet())
            {
                GM_LOGGER.debug("soft filter {}: count({})", entry.getKey(), entry.getValue());
            }
        }
    }

    private void processVariant(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        GM_LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

        ++mProcessedVariants;

        if(mProcessedVariants > 0 && (mProcessedVariants % 100000) == 0)
        {
            GM_LOGGER.debug("sample({}) processed {} variants, VCF-unmatched({})",
                    mConfig.SampleId, mProcessedVariants, mSvFactory.unmatched().size());
        }

        // first check hard-filters which only operate on raw variant context info
        String mateId = StructuralVariantFactory.mateId(variant);

        if(mateId != null)
        {
            if(mHardFilteredVcfIds.contains(mateId))
                return; // already filtered
        }

        if(mHardFilters.isFiltered(variant, genotypeIds))
        {
            if(mateId != null)
            {
                // other breakend not already filtered otherwise will have exited earlier, so either cached or not seen yet
                // no point in keeping the other breakend if cached
                if(mSvFactory.hasUnmatchedVariant(mateId))
                {
                    mSvFactory.removeUnmatchedVariant(mateId);
                    return;
                }

                mHardFilteredVcfIds.add(variant.getID());
            }
            else
            {
                // mate ID not set or a single breakend, but either way no need to register hard-filtered ID
            }

            registerFilter(FilterType.HARD_FILTERED);
            return;
        }

        int currentSvCount = mSvFactory.results().size();
        mSvFactory.addVariantContext(variant);

        // wait for both breakends to be added
        if(currentSvCount == mSvFactory.results().size())
            return;

        final StructuralVariant sv = popLastSv(); // get and clear from storage

        if(sv == null)
            return;

        // check hard filters again on the pair of breakends to ensure they both match a known pair
        if(!mHardFilters.keepHotspotVariant(sv))
        {
            registerFilter(FilterType.HARD_FILTERED);
            return;
        }

        // optionally filter out by config
        if(mConfig.excludeVariant(sv))
            return;

        SvData svData = new SvData(sv, genotypeIds);
        mSvDataCache.addSvData(svData);
    }

    private void registerFilter(final FilterType type)
    {
        Integer count = mFilterCounts.get(type);
        mFilterCounts.put(type, count != null ? count + 1 : 1);
    }

    private final StructuralVariant popLastSv()
    {
        if(mSvFactory.results().isEmpty())
            return null;

        StructuralVariant sv = mSvFactory.results().get(0);
        mSvFactory.results().remove(0);

        return sv;
    }

    private void clearSampleData()
    {
        /*
        mProcessedVariants = 0;
        mSampleSvData.clear();
        mSvFactory = new StructuralVariantFactory(new AlwaysPassFilter());
        mLinkAnalyser.clear();
         */
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        GripssConfig.addCommandLineOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        GripssApplication gripssApplication = new GripssApplication(cmd);
        gripssApplication.run();

        GM_LOGGER.info("VCF processing complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
