package com.hartwig.hmftools.amber.utils;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class Mappability
{
    private BufferedReader mFileReader;
    private String mCurrentChromosome;
    private List<MapEntry> mChromosomeEntries;
    private MapEntry mNextChromosomeEntry;
    private int mCurrentIndex;
    private boolean mHasValidData;

    public static final double INVALID = 0;

    public static final String MAPPABILITY_BED = "mappability_bed";

    public Mappability(final ConfigBuilder configBuilder)
    {
        mFileReader = null;
        mCurrentChromosome = "";
        mChromosomeEntries = Lists.newArrayList();
        mNextChromosomeEntry = null;
        mCurrentIndex = 0;
        mHasValidData = true;

        if(configBuilder.hasValue(MAPPABILITY_BED))
        {
            initialiseFile(configBuilder.getValue(MAPPABILITY_BED));
        }
    }

    public boolean hasData() { return mFileReader != null; }
    public boolean hasValidData() { return mHasValidData; }

    public double getMappability(final String chromosome, final int position)
    {
        if(!mHasValidData || (mFileReader == null && mChromosomeEntries.isEmpty()))
            return INVALID;

        while(!mChromosomeEntries.isEmpty() || mFileReader != null)
        {
            if(chromosome.equals(mCurrentChromosome))
                break;

            loadEntries(chromosome);
        }

        return findMappability(chromosome, position);
    }

    private double findMappability(final String chromosome, final int position)
    {
        if(!chromosome.equals(mCurrentChromosome))
            return INVALID;

        for(int i = mCurrentIndex; i < mChromosomeEntries.size(); ++i)
        {
            MapEntry entry = mChromosomeEntries.get(i);

            if(entry.Region.containsPosition(position))
            {
                mCurrentIndex = i;
                return entry.Mappability;
            }

            // take previous if the next is past this variant
            if(position < entry.Region.start() && i > 0)
            {
                MapEntry prevEntry = mChromosomeEntries.get(i - 1);
                return prevEntry.Mappability;
            }
        }

        return INVALID;
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
        }
        catch(IOException e)
        {
            AMB_LOGGER.error("failed to load mappability file({}): {}", filename, e.toString());
            mHasValidData = false;
        }
    }

    private void loadEntries(final String requestedChromosome)
    {
        if(mFileReader == null)
        {
            mChromosomeEntries.clear();
            return;
        }

        // clear all but the last if it's a match for the required chromosome
        mChromosomeEntries.clear();

        if(mNextChromosomeEntry != null && mNextChromosomeEntry.Region.chromosome().equals(requestedChromosome))
            mChromosomeEntries.add(mNextChromosomeEntry);

        mNextChromosomeEntry = null;

        mCurrentChromosome = requestedChromosome;
        mCurrentIndex = 0;

        try
        {
            String line = null;
            boolean foundRequested = !mChromosomeEntries.isEmpty();

            while((line = mFileReader.readLine()) != null)
            {
                final String[] values = line.split(TSV_DELIM, -1); // eg: 1       0       10000   0.000000

                if(values.length != 4)
                {
                    mHasValidData = false;
                    mFileReader = null;
                    return;
                }

                String chromosome = values[0];

                if(chromosome.equals(requestedChromosome))
                {
                    foundRequested = true;
                }
                else
                {
                    if(!foundRequested)
                        continue; // skip past this until a match is found
                }

                MapEntry entry = new MapEntry(
                        new ChrBaseRegion(chromosome, Integer.parseInt(values[1]) + 1, Integer.parseInt(values[2])),
                        Double.parseDouble(values[3]));

                if(chromosome.equals(requestedChromosome))
                {
                    mChromosomeEntries.add(entry);
                }
                else
                {
                    mNextChromosomeEntry = entry;
                    break;
                }
            }

            if(line == null)
                mFileReader = null;
        }
        catch(IOException e)
        {
            AMB_LOGGER.error("failed to load mappability file: {}", e.toString());
            mHasValidData = false;
        }
    }

    private class MapEntry
    {
        public final ChrBaseRegion Region;
        public final double Mappability;

        public MapEntry(final ChrBaseRegion region, final double mappability)
        {
            Region = region;
            Mappability = mappability;
        }

        public String toString() { return String.format("%s map(%.4f)", Region, Mappability); }
    }
}
