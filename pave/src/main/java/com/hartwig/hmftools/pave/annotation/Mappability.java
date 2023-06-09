package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.pave.VariantData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class Mappability
{
    private BufferedReader mFileReader;
    private String mCurrentChromosome;
    private int mCurrentIndex;
    private List<MapEntry> mEntries;
    private MapEntry mNextChromosomeEntry;
    private boolean mHasValidData;

    private static final int MAX_LIST_COUNT = 1000;
    public static final String MAPPABILITY_BED = "mappability_bed";

    public static final String MAPPABILITY = "MAPPABILITY";
    public static final String MAPPABILITY_DESC = "GEM mappability in 150 base window";

    public Mappability(final CommandLine cmd)
    {
        mFileReader = null;
        mCurrentChromosome = "";
        mCurrentIndex = 0;
        mEntries = Lists.newArrayList();
        mNextChromosomeEntry = null;
        mHasValidData = true;

        if(cmd.hasOption(MAPPABILITY_BED))
        {
            initialiseFile(cmd.getOptionValue(MAPPABILITY_BED));
        }
    }

    public boolean hasData() { return mFileReader != null; }
    public boolean hasValidData() { return mHasValidData; }

    public void annotateVariant(final VariantData variant)
    {
        if(mFileReader == null && mEntries.isEmpty())
            return;

        while(!mEntries.isEmpty() || mFileReader != null)
        {
            if(checkEntries(variant))
                return;

            loadEntries(variant.Chromosome);
        }

        PV_LOGGER.warn("variant({}) no mappability entry found", variant);
    }

    private boolean checkEntries(final VariantData variant)
    {
        if(!variant.Chromosome.equals(mCurrentChromosome))
            return false;

        if(mEntries.isEmpty() || mEntries.get(mEntries.size() - 1).Region.end() < variant.Position)
            return false;

        for(int i = mCurrentIndex; i < mEntries.size(); ++i)
        {
            MapEntry entry = mEntries.get(i);

            if(entry.Region.containsPosition(variant.Position))
            {
                setMappability(variant, entry.Mappability);
                mCurrentIndex = i;
                return true;
            }

            // take previous if the next is past this variant
            if(variant.Position < entry.Region.start() && i > 0)
            {
                MapEntry prevEntry = mEntries.get(i - 1);
                setMappability(variant, prevEntry.Mappability);
                return true;
            }
        }

        return false;
    }

    private void setMappability(final VariantData variant, double mappability)
    {
        if(!variant.context().getCommonInfo().hasAttribute(MAPPABILITY))
        {
            variant.context().getCommonInfo().putAttribute(MAPPABILITY, mappability);
        }
    }

    public static void addHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(MAPPABILITY, 1, VCFHeaderLineType.Float, MAPPABILITY_DESC));
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(MAPPABILITY_BED, true, "Mappability BED file");
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
                PV_LOGGER.info("loaded mappability file({})", filename);
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
        {
            mEntries.clear();
            return;
        }

        // clear all but the last if it's a match for the required chromosome
        mEntries.clear();

        if(mNextChromosomeEntry != null && mNextChromosomeEntry.Region.chromosome().equals(requestedChromosome))
            mEntries.add(mNextChromosomeEntry);

        mNextChromosomeEntry = null;

        mCurrentChromosome = requestedChromosome;
        mCurrentIndex = 0;

        try
        {
            String line = null;
            boolean foundRequested = !mEntries.isEmpty();

            while((line = mFileReader.readLine()) != null)
            {
                final String[] values = line.split("\t", -1); // eg: 1       0       10000   0.000000

                if(values.length != 4)
                {
                    mHasValidData = false;
                    mFileReader = null;

                    PV_LOGGER.error("invalid mappability entry({})", line);
                    System.exit(1);
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
                    mEntries.add(entry);

                    if(mEntries.size() >= MAX_LIST_COUNT)
                        return;
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
            PV_LOGGER.error("failed to load mappability file: {}", e.toString());
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
