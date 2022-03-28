package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.ITEM_DELIM;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.VariantTier;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class PonAnnotation
{
    private final Map<String,Map<Integer,List<PonVariantData>>> mPonEntries; // mapped by chromosome then position
    private final boolean mLoadOnDemand;

    private final String mPonFilename;
    private BufferedReader mFileReader;
    private String mCurrentChromosome;
    private int mColumnCount;
    private boolean mHasValidData;

    private final Map<VariantTier,PonFilters> mPonFilters;

    public static final String PON_COUNT = "PON_COUNT";
    public static final String PON_MAX = "PON_MAX";
    public static final String PON_FILTER = "PON";
    public static final String PON_ARTEFACT_FILTER = "PONArtefact";

    public PonAnnotation(final String filename, boolean loadOnDemand)
    {
        mPonEntries = Maps.newHashMap();
        mLoadOnDemand = loadOnDemand;
        mPonFilename = filename;
        mFileReader = null;
        mCurrentChromosome = "";
        mColumnCount = -1;
        mHasValidData = false;

        if(filename != null && !filename.isEmpty())
        {
            loadPonFile(filename);
        }

        mPonFilters = Maps.newHashMap();
    }

    public boolean hasValidData() { return mHasValidData; }

    public boolean isEnabled() { return mPonFilename != null; }

    public boolean loadFilters(final String filtersConfig)
    {
        if(filtersConfig == null)
            return true;

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

    public void annotateVariant(final VariantData variant)
    {
        PonVariantData ponData = getPonData(variant);
        if(ponData == null)
            return;

        variant.setPonFrequency(ponData.Samples, ponData.MaxSampleReads);

        VariantTier tier = variant.tier();

        PonFilters filters = mPonFilters.containsKey(tier) ? mPonFilters.get(tier) : mPonFilters.get(VariantTier.UNKNOWN);

        if(filters == null)
            return;

        if(ponData.Samples >= filters.RequiredSampleCount && ponData.MaxSampleReads >= filters.RequiredMaxReadCount)
            variant.addFilter(PON_FILTER);
    }

    public PonVariantData getPonData(final VariantData variant)
    {
        if(mLoadOnDemand && !variant.Chromosome.equals(mCurrentChromosome))
        {
            mPonEntries.remove(mCurrentChromosome);
            mCurrentChromosome = variant.Chromosome;
            loadPonEntries();
        }
        Map<Integer,List<PonVariantData>> posMap = mPonEntries.get(variant.Chromosome);

        if(posMap == null)
            return null;

        List<PonVariantData> posList = posMap.get(variant.Position);

        if(posList == null)
            return null;

        return posList.stream().filter(x -> x.matches(variant.Ref, variant.Alt)).findFirst().orElse(null);
    }

    private void loadPonFile(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            mFileReader = createBufferedReader(filename);

            String line = mFileReader.readLine();
            final String[] values = line.split("\t", -1);
            mColumnCount = values.length;

            if(mColumnCount < 4)
            {
                PV_LOGGER.error("pon file({}) has insufficient column count({})", filename, mColumnCount);
                mFileReader = null;
            }

            if(!mLoadOnDemand)
            {
                // load the full file
                loadPonEntries();
            }

            mHasValidData = true;
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load PON file({}): {}", filename, e.toString());
            mHasValidData = false;
        }
    }

    private void loadPonEntries()
    {
        if(mFileReader == null)
            return;

        try
        {
            int currentPos = 0;

            String requestedChromosome = mLoadOnDemand ? mCurrentChromosome : null;

            Map<Integer,List<PonVariantData>> posMap = null;
            List<PonVariantData> posList = null;

            String line = null;
            String currentChr = "";
            boolean foundRequested = false;

            if(requestedChromosome != null && mPonEntries.containsKey(requestedChromosome))
            {
                currentChr = requestedChromosome;
                posMap = mPonEntries.get(requestedChromosome);
            }

            int itemCount = 0;

            while((line = mFileReader.readLine()) != null)
            {
                final String[] values = line.split("\t", -1);

                int colIndex = 0;
                String chromosome = values[colIndex++];

                if(requestedChromosome != null)
                {
                    if(chromosome.equals(requestedChromosome))
                    {
                        foundRequested = true;
                    }
                    else
                    {
                        if(!foundRequested)
                            continue; // skip past this until a match is found

                        // otherwise record this entry and then break
                    }
                }

                int position = Integer.parseInt(values[colIndex++]);
                String ref = values[colIndex++];
                String alt = values[colIndex++];

                int sampleCount = mColumnCount > colIndex ? Integer.parseInt(values[colIndex++]) : 0;
                int maxReadsCount = mColumnCount > colIndex ? Integer.parseInt(values[colIndex++]) : 0;
                int totalReadsCount = mColumnCount > colIndex ? Integer.parseInt(values[colIndex++]) : 0;

                if(!chromosome.equals(currentChr))
                {
                    currentChr = chromosome;
                    currentPos = position;
                    posMap = Maps.newHashMap();
                    mPonEntries.put(chromosome, posMap);

                    posList = Lists.newArrayList();
                    posMap.put(position, posList);
                }
                else if(currentPos != position)
                {
                    currentPos = position;
                    posList = Lists.newArrayList();
                    posMap.put(position, posList);
                }

                posList.add(new PonVariantData(ref, alt, sampleCount, maxReadsCount, totalReadsCount));
                ++itemCount;

                if(requestedChromosome != null && foundRequested && !chromosome.equals(requestedChromosome))
                {
                    --itemCount;
                    break;
                }
            }

            if(requestedChromosome != null)
            {
                PV_LOGGER.debug("pon file({}) loaded {} entries for chromosome({})", mPonFilename, itemCount, requestedChromosome);
            }
            else
            {
                PV_LOGGER.info("pon file({}) loaded {} entries", mPonFilename, itemCount);
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load PON file: {}", e.toString());
            mHasValidData = false;
            return;
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

    private class PonFilters
    {
        public final int RequiredSampleCount;
        public final int RequiredMaxReadCount;

        public PonFilters(final int requiredSampleCount, final int requiredMaxReadCount)
        {
            RequiredSampleCount = requiredSampleCount;
            RequiredMaxReadCount = requiredMaxReadCount;
        }
    }

}
