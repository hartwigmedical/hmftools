package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.MAPPABILITY;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.MAPPABILITY_DESC;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class Mappability extends AnnotationData implements Callable<Void>
{
    private BufferedReader mFileReader;
    private final Map<String, MappabilityChrCache> mChrCacheMap;
    private boolean mHasValidData;

    public static final String MAPPABILITY_BED = "mappability_bed";

    public Mappability(final ConfigBuilder configBuilder)
    {
        mFileReader = null;
        mHasValidData = true;
        mChrCacheMap = Maps.newHashMap();

        if(configBuilder.hasValue(MAPPABILITY_BED))
        {
            initialiseFile(configBuilder.getValue(MAPPABILITY_BED));
        }
    }

    @Override
    public String type() { return "Mappability"; }

    @Override
    public boolean enabled() { return mFileReader != null; }

    @Override
    public boolean hasValidData() { return mHasValidData; }

    public synchronized MappabilityChrCache getChromosomeCache(final String chromosome)
    {
        MappabilityChrCache chrCache = mChrCacheMap.get(chromosome);

        if(chrCache != null && chrCache.isComplete())
            return chrCache;

        loadEntries(chromosome);
        return mChrCacheMap.get(chromosome);
    }

    @Override
    public synchronized void onChromosomeComplete(final String chromosome)
    {
        MappabilityChrCache chrCache = mChrCacheMap.get(chromosome);

        if(chrCache != null)
        {
            chrCache.clear();
            mChrCacheMap.remove(chromosome);
        }
    }

    @Override
    public Void call()
    {
        for(String chromosome : mInitialChromosomes)
        {
            loadEntries(chromosome);
        }

        return null;
    }

    public static void addHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(MAPPABILITY, 1, VCFHeaderLineType.Float, MAPPABILITY_DESC));
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(MAPPABILITY_BED, false, "Mappability BED file");
    }

    private void initialiseFile(final String filename)
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

            if(mFileReader != null)
                PV_LOGGER.debug("loading mappability file({})", filename);
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load mappability file({}): {}", filename, e.toString());
            mHasValidData = false;
        }
    }

    private void loadEntries(final String requestedChromosome)
    {
        if(mFileReader == null)
            return;

        try
        {
            String line = null;
            MappabilityChrCache currentCache = mChrCacheMap.get(requestedChromosome);
            boolean foundRequested = currentCache != null;

            while((line = mFileReader.readLine()) != null)
            {
                final String[] values = line.split(TSV_DELIM, -1); // eg: 1       0       10000   0.000000

                if(values.length != 4)
                {
                    mHasValidData = false;
                    mFileReader = null;

                    PV_LOGGER.error("invalid mappability entry({})", line);
                    System.exit(1);
                }

                String chromosome = values[0];

                boolean exitOnNew = false;

                if(currentCache == null || !currentCache.Chromosome.equals(chromosome))
                {
                    if(currentCache != null)
                    {
                        PV_LOGGER.debug("chr({}) loaded {} mappability entries", currentCache.Chromosome, currentCache.entryCount());
                        currentCache.setComplete();
                    }

                    currentCache = mChrCacheMap.get(chromosome);

                    if(currentCache == null)
                    {
                        currentCache = new MappabilityChrCache(chromosome);
                        mChrCacheMap.put(chromosome, currentCache);
                    }

                    if(!foundRequested)
                    {
                        foundRequested = chromosome.equals(requestedChromosome);
                    }
                    else if(HumanChromosome.chromosomeRank(chromosome) > HumanChromosome.chromosomeRank(requestedChromosome))
                    {
                        exitOnNew = true;
                    }
                }

                currentCache.addEntry(Integer.parseInt(values[1]) + 1, Integer.parseInt(values[2]), Double.parseDouble(values[3]));

                if(exitOnNew)
                    break;
            }

            if(line == null)
                mFileReader = null;
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load mappability file: {}", e.toString());
            mHasValidData = false;
        }
    }
}
