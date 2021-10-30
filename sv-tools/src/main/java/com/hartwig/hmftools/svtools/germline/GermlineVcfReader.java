package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.svtools.germline.GermlineUtils.GM_LOGGER;
import static com.hartwig.hmftools.svtools.germline.GermlineUtils.findVcfFiles;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.loadVcfFiles;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;

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

public class GermlineVcfReader
{
    private final GermlineFilters mFilter;
    private final GermlineVcfConfig mConfig;
    public final FilterConstants mFilterConstants;

    private final LinkAnalyser mLinkAnalyser;
    private final PonCache mPonCache;
    private final HotspotCache mHotspotCache;
    private final HardFilters mHardFilters;

    private int mProcessedVariants;
    private StructuralVariantFactory mSvFactory;
    private final Set<String> mHardFilteredVcfIds;
    private final List<SvData> mSampleSvData;
    private final Map<FilterType,Integer> mFilterCounts;

    public GermlineVcfReader(final CommandLine cmd)
    {
        mConfig = new GermlineVcfConfig(cmd);
        mFilterConstants = FilterConstants.from(cmd);

        mFilter = new GermlineFilters(mConfig);

        mSvFactory = new StructuralVariantFactory(new AlwaysPassFilter());
        mHardFilteredVcfIds = Sets.newHashSet();

        mLinkAnalyser = new LinkAnalyser();

        GM_LOGGER.info("loading reference data");
        mPonCache = new PonCache(cmd, mConfig.RestrictedChromosomes);
        mHotspotCache = new HotspotCache(cmd);

        mHardFilters = new HardFilters(mFilterConstants, mHotspotCache);

        mProcessedVariants = 0;
        mSampleSvData = Lists.newArrayList();
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

            clearSampleData();

            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    vcfFile, new VCFCodec(), false);

            GenotypeIds genotypeIds = VcfUtils.parseVcfSampleIds((VCFHeader)reader.getHeader(), mConfig.ReferenceId, mConfig.SampleId);

            if(genotypeIds == null)
            {
                System.exit(1);
            }

            reader.iterator().forEach(x -> processVariant(x, genotypeIds));

            GM_LOGGER.info("sample({}) read VCF: unmatched({}) results({})",
                    mConfig.SampleId, mSvFactory.unmatched().size(), mSvFactory.results().size());

            if(mSampleSvData.isEmpty())
                return;

            for(final SvData svData : mSampleSvData)
            {
                // annotate with the PON
                svData.setPonCount(mPonCache.getPonCount(svData));

                //mLinkAnalyser.cacheAssemblyData(svData);
                //mLinkAnalyser.populateAssemblyLinks(svData);
                // writeCsv(germlineSV);
            }
        }
        catch(IOException e)
        {
            GM_LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }

        int hardFiltered = mFilterCounts.values().stream().mapToInt(x -> x.intValue()).sum();

        GM_LOGGER.info("sample({}) read {} variants: hardFiltered({}) cached({})",
                mConfig.SampleId, mProcessedVariants, hardFiltered, mSampleSvData.size());
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

        if(mHardFilters.failsMinTumorQuality(variant, genotypeIds.TumorOrdinal))
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

            registerFilter(FilterType.HARD_MIN_QUAL);
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
        FilterType hardFilter = mHardFilters.getFilterType(sv);

        if(hardFilter != FilterType.PASS)
        {
            registerFilter(hardFilter);
            return;
        }

        // optionally filter out by config
        if(mConfig.excludeVariant(sv))
            return;

        SvData svData = SvData.from(sv, genotypeIds);
        mSampleSvData.add(svData);
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
        GermlineVcfConfig.addCommandLineOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        GermlineVcfReader germlineVcfReader = new GermlineVcfReader(cmd);
        germlineVcfReader.run();

        GM_LOGGER.info("VCF processing complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
