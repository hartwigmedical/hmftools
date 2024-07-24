package com.hartwig.hmftools.pave.annotation;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.PON_FILTER_HOTSPOT_MAX_READS;
import static com.hartwig.hmftools.pave.PaveConstants.PON_FILTER_HOTSPOT_SAMPLE_COUNT;
import static com.hartwig.hmftools.pave.PaveConstants.PON_FILTER_OTHER_TIER_MAX_READS;
import static com.hartwig.hmftools.pave.PaveConstants.PON_FILTER_OTHER_TIER_SAMPLE_COUNT;
import static com.hartwig.hmftools.pave.PaveConstants.PON_FILTER_PANEL_MAX_READS;
import static com.hartwig.hmftools.pave.PaveConstants.PON_FILTER_PANEL_SAMPLE_COUNT;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.StringCache;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.pave.VariantData;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class PonAnnotation extends AnnotationData implements Callable
{
    private final String mPonFilename;
    private BufferedReader mFileReader;
    private int mColumnCount;
    private boolean mHasValidData;

    private final Map<VariantTier,PonFilters> mPonFilters;
    private final Map<String,PonChrCache> mChrCacheMap;
    private final StringCache mStringCache;

    public static final String PON_COUNT = "PON_COUNT";
    public static final String PON_MAX = "PON_MAX";
    public static final String PON_FILTER = "PON";
    public static final String PON_ARTEFACT_FILTER = "PONArtefact";

    public PonAnnotation(final String filename, boolean loadOnDemand)
    {
        mPonFilename = filename;
        mFileReader = null;
        mColumnCount = -1;
        mHasValidData = true;
        mChrCacheMap = Maps.newHashMap();
        mStringCache = new StringCache();

        if(filename != null && !filename.isEmpty())
        {
            initialiseFile(filename, loadOnDemand);
        }

        mPonFilters = Maps.newHashMap();

        // set defaults
        mPonFilters.put(VariantTier.HOTSPOT, new PonFilters(PON_FILTER_HOTSPOT_SAMPLE_COUNT, PON_FILTER_HOTSPOT_MAX_READS));
        mPonFilters.put(VariantTier.PANEL, new PonFilters(PON_FILTER_PANEL_SAMPLE_COUNT, PON_FILTER_PANEL_MAX_READS));
        mPonFilters.put(VariantTier.UNKNOWN, new PonFilters(PON_FILTER_OTHER_TIER_SAMPLE_COUNT, PON_FILTER_OTHER_TIER_MAX_READS));
    }

    @Override
    public String type() { return "PON"; }

    @Override
    public boolean enabled() { return mPonFilename != null; }

    @Override
    public boolean hasValidData() { return mHasValidData; }

    public boolean loadFilters(final String filtersConfig)
    {
        if(filtersConfig == null)
            return true;

        mPonFilters.clear();

        String[] tierFilters = filtersConfig.split(ITEM_DELIM, -1);

        try
        {
            for(String tierFilter : tierFilters)
            {
                String[] values = tierFilter.split(":", -1);
                VariantTier tier = VariantTier.fromString(values[0]);
                int reqSampleCount = Integer.parseInt(values[1]);
                int reqMaxReads = Integer.parseInt(values[2]);
                mPonFilters.put(tier, new PonFilters(reqSampleCount, reqMaxReads));

                PV_LOGGER.info("loaded PON filter: tier({}) sampleCount({}) maxReads({})",
                        tier, reqSampleCount, reqMaxReads);
            }
        }
        catch(Exception e)
        {
            PV_LOGGER.error("invalid PON filters(): {}", filtersConfig, e.toString());
            mHasValidData = false;
            return false;
        }

        return true;
    }

    public synchronized PonChrCache getChromosomeCache(final String chromosome)
    {
        PonChrCache chrCache = mChrCacheMap.get(chromosome);

        if(chrCache != null && chrCache.isComplete())
            return chrCache;

        loadPonEntries(chromosome);
        return mChrCacheMap.get(chromosome);
    }

    @Override
    public synchronized void onChromosomeComplete(final String chromosome)
    {
        PonChrCache chrCache = mChrCacheMap.get(chromosome);

        if(chrCache != null)
        {
            chrCache.clear();
            mChrCacheMap.remove(chromosome);
        }
    }

    @Override
    public Long call()
    {
        for(String chromosome : mInitialChromosomes)
        {
            loadPonEntries(chromosome);
        }

        return (long)0;
    }

    public void annotateVariant(final VariantData variant, final PonChrCache chrCache)
    {
        PonVariantData ponData = chrCache.getPonData(variant);
        if(ponData == null)
            return;

        variant.setPonFrequency(ponData.Samples, ponData.MaxSampleReads, ponData.meanReadCount());
    }

    public boolean filterOnTierCriteria(final VariantTier tier, final int ponSampleCount, final int ponMaxSampleReads)
    {
        PonFilters filters = mPonFilters.containsKey(tier) ? mPonFilters.get(tier) : mPonFilters.get(VariantTier.UNKNOWN);

        if(filters == null)
            return false;

        return ponSampleCount >= filters.SampleCount && ponMaxSampleReads >= filters.MaxReadCount;
    }

    public boolean hasEntry(final String chromosome, final int position, final String ref, final String alt)
    {
        PonChrCache chrCache = mChrCacheMap.get(chromosome);
        return chrCache != null ? chrCache.hasEntry(position, ref, alt) : false;
    }

    private void initialiseFile(final String filename, boolean loadOnDemand)
    {
        if(filename == null)
            return;

        if(!Files.exists(Paths.get(filename)))
        {
            mHasValidData = false;
            return;
        }

        try
        {
            mFileReader = createBufferedReader(filename);

            String line = mFileReader.readLine();
            final String[] values = line.split(TSV_DELIM, -1);
            mColumnCount = values.length;

            if(mColumnCount < 4)
            {
                PV_LOGGER.error("pon file({}) has insufficient column count({})", filename, mColumnCount);
                mFileReader = null;
            }

            mHasValidData = true;

            if(!loadOnDemand)
                loadPonEntries(null);

        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load PON file({}): {}", filename, e.toString());
            mHasValidData = false;
        }
    }

    private void loadPonEntries(@Nullable final String requestedChromosome)
    {
        if(mFileReader == null)
            return;

        try
        {
            String line = null;
            PonChrCache currentCache = requestedChromosome != null ? mChrCacheMap.get(requestedChromosome) : null;
            boolean foundRequested = currentCache != null;
            int itemCount = 0;

            while((line = mFileReader.readLine()) != null)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                int colIndex = 0;
                String chromosome = values[colIndex++];

                boolean exitOnNew = false;

                if(currentCache == null || !currentCache.Chromosome.equals(chromosome))
                {
                    if(currentCache != null && !currentCache.isComplete())
                    {
                        PV_LOGGER.debug("chr({}) loaded {} PON entries", currentCache.Chromosome, currentCache.entryCount());

                        currentCache.setComplete();
                    }

                    currentCache = mChrCacheMap.get(chromosome);

                    if(currentCache == null)
                    {
                        currentCache = new PonChrCache(chromosome, mStringCache);
                        mChrCacheMap.put(chromosome, currentCache);
                    }

                    if(requestedChromosome != null)
                    {
                        if(!foundRequested)
                        {
                            foundRequested = chromosome.equals(requestedChromosome);
                        }
                        else if(HumanChromosome.chromosomeRank(chromosome) > HumanChromosome.chromosomeRank(requestedChromosome))
                        {
                            exitOnNew = true;
                        }
                    }
                }

                int position = Integer.parseInt(values[colIndex++]);
                String ref = values[colIndex++];
                String alt = values[colIndex++];

                int sampleCount = mColumnCount > colIndex ? Integer.parseInt(values[colIndex++]) : 0;
                int maxReadsCount = mColumnCount > colIndex ? Integer.parseInt(values[colIndex++]) : 0;
                int totalReadsCount = mColumnCount > colIndex ? Integer.parseInt(values[colIndex++]) : 0;

                currentCache.addEntry(position, ref, alt, sampleCount, maxReadsCount, totalReadsCount);
                ++itemCount;

                if(exitOnNew)
                    break;
            }

            if(requestedChromosome == null)
            {
                PV_LOGGER.info("pon file({}) loaded {} entries", mPonFilename, itemCount);
            }

            if(line == null)
                mFileReader = null;
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load PON file: {}", e.toString());
            mHasValidData = false;
        }
    }

    public static void addHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(
                PON_COUNT, 1, VCFHeaderLineType.Integer, "Cohort frequency for variant"));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                PON_MAX, 1, VCFHeaderLineType.Integer, "Max read depth in any sample with variant"));

        header.addMetaDataLine(new VCFFilterHeaderLine(PON_ARTEFACT_FILTER, "Filter PON artefact"));
        header.addMetaDataLine(new VCFFilterHeaderLine(PON_FILTER, "Filter PON variant"));
    }

    private static class PonFilters
    {
        public final int SampleCount;
        public final int MaxReadCount;

        public PonFilters(final int sampleCount, final int maxReadCount)
        {
            SampleCount = sampleCount;
            MaxReadCount = maxReadCount;
        }

        public String toString() { return format("samples(%d) maxReads(%d)", SampleCount, MaxReadCount); }
    }
}
