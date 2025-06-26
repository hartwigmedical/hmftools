package com.hartwig.hmftools.common.variant.pon;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ_DESC;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.perf.StringCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GnomadCache
{
    private final Map<String,GnomadChrCache> mChrCacheMap;
    private final RefGenomeVersion mRefGenomeVersion;
    private final Map<String,String> mChromosomeFiles;
    private boolean mHasValidData;
    private final boolean mEnabled;
    private final String mGnomadFilename;
    private final StringCache mStringCache;

    public static final String GNOMAD_FREQUENCY_FILE = "gnomad_freq_file";
    public static final String GNOMAD_FREQUENCY_DIR = "gnomad_freq_dir";

    public static final String PON_GNOMAD_FILTER = "PONGnomad";

    public static final String GNOMAD_FILE_ID = "gnomad_variants";

    protected static final Logger LOGGER = LogManager.getLogger(GnomadCache.class);

    public GnomadCache(final RefGenomeVersion refGenomeVersion, @Nullable final String gnomadFilename, @Nullable final String gnomadDirectory)
    {
        mChrCacheMap = Maps.newHashMap();
        mChromosomeFiles = Maps.newHashMap();
        mHasValidData = true;
        mStringCache = new StringCache();

        mRefGenomeVersion = refGenomeVersion;

        if(gnomadFilename != null)
        {
            mEnabled = true;
            mGnomadFilename = gnomadFilename;
        }
        else if(gnomadDirectory != null)
        {
            mEnabled = true;
            mGnomadFilename = null;
            String gnomadDir = gnomadDirectory;
            mHasValidData = loadAllFrequencyFiles(gnomadDir);
        }
        else
        {
            mGnomadFilename = null;
            mEnabled = false;
        }
    }

    public boolean enabled() { return mGnomadFilename != null || !mChromosomeFiles.isEmpty(); }
    public boolean hasValidData() { return mHasValidData; }

    public void initialise(final List<String> initialChromosomes)
    {
        if(mGnomadFilename != null)
        {
            loadChromosomeEntries(mGnomadFilename, null);
        }
        else if(!mChromosomeFiles.isEmpty())
        {
            for(String chromosome : initialChromosomes)
            {
                getChromosomeCache(chromosome);
            }
        }
    }

    public synchronized GnomadChrCache getChromosomeCache(final String chromosome)
    {
        if(!mEnabled)
            return null;

        GnomadChrCache chrCache = mChrCacheMap.get(chromosome);

        if(chrCache != null)
            return chrCache;

        String chrFilename = mChromosomeFiles.get(chromosome);

        if(chrFilename == null)
        {
            LOGGER.warn("missing Gnomad file for chromosome({})", chromosome);
            return null;
        }

        loadChromosomeEntries(chrFilename, chromosome);
        return mChrCacheMap.get(chromosome);
    }

    public synchronized void removeCompleteChromosome(final String chromosome)
    {
        GnomadChrCache chrCache = mChrCacheMap.get(chromosome);

        if(chrCache != null)
        {
            chrCache.clear();
            mChrCacheMap.remove(chromosome);
        }
    }

    public static String formFileId(final String dir, final String chromosome, final String outputId)
    {
        String outputFile = dir + GNOMAD_FILE_ID;

        if(chromosome != null && !chromosome.isEmpty())
            outputFile += "_chr" + chromosome;

        if(outputId != null)
            outputFile += "_" + outputId;

        outputFile += ".csv";
        return outputFile;
    }

    private boolean loadChromosomeEntries(final String filename, final String fileChromosome)
    {
        // if file chromosome is supplied then it is not read from the input file
        if(filename == null)
            return false;

        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            int itemCount = 0;
            String line = fileReader.readLine(); // skip header
            String currentChr = "";
            GnomadChrCache currentChrCache = null;

            Integer chrIndex = fileChromosome != null ? -1 : 0;
            int index = fileChromosome != null ? 0 : 1;
            int posIndex = index++;
            int refIndex = index++;
            int altIndex = index++;
            int freqIndex = index;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(CSV_DELIM, -1);

                String chromosome = fileChromosome != null ? fileChromosome : values[chrIndex];
                int position = Integer.parseInt(values[posIndex]);
                String ref = values[refIndex];
                String alt = values[altIndex];
                double frequency = Double.parseDouble(values[freqIndex]);

                if(!chromosome.equals(currentChr))
                {
                    currentChr = chromosome;
                    currentChrCache = new GnomadChrCache(chromosome, mStringCache);
                    mChrCacheMap.put(chromosome, currentChrCache);
                }

                currentChrCache.addEntry(position, ref, alt,frequency);

                ++itemCount;
            }

            if(mChromosomeFiles.isEmpty())
            {
                LOGGER.info("loaded {} Gnomad frequency records from file({})", itemCount, filename);
            }
            else if(currentChrCache != null)
            {
                LOGGER.debug("chr({}) loaded {} Gnomad frequency records",
                        currentChrCache.Chromosome, currentChrCache.entryCount());
            }

            return true;
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load Gnomad frequency file({}): {}", filename, e.toString());
            return false;
        }
    }

    private boolean loadAllFrequencyFiles(final String gnomadDir)
    {
        try
        {
            final Stream<Path> stream = Files.walk(Paths.get(gnomadDir), 1, FileVisitOption.FOLLOW_LINKS);

            List<String> files = stream.filter(x -> !x.toFile().isDirectory())
                    .map(x -> x.toFile().toString())
                    .filter(x -> x.contains(GNOMAD_FILE_ID))
                    .collect(Collectors.toList());

            for(HumanChromosome humanChr : HumanChromosome.values())
            {
                String fileChrStrNoId = formFileId(gnomadDir, humanChr.toString(), null);
                String fileChrStrWithId = GNOMAD_FILE_ID + "_chr" + humanChr + "_";

                // expect file name: gnomad_variants_chr10_v38.csv.gz

                String chrFile = files.stream()
                        .filter(x -> x.endsWith(fileChrStrNoId) || x.contains(fileChrStrWithId))
                        .findFirst().orElse(null);

                String chrStr = mRefGenomeVersion.versionedChromosome(humanChr.toString());

                if(chrFile == null)
                {
                    LOGGER.error("missing Gnomad chromosome({}) file", chrStr);
                    return false;
                }

                mChromosomeFiles.put(chrStr, chrFile);
            }

            return true;
        }
        catch(IOException e)
        {
            LOGGER.error("failed to find Gnomad chromosome files in dir({}): {}", gnomadDir, e.toString());
            return false;
        }
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(GNOMAD_FREQUENCY_FILE, false, "Gnomad frequency file");
        configBuilder.addPath(GNOMAD_FREQUENCY_DIR, false, "Gnomad frequency directory");
    }

    public static void addAnnotationHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(GNOMAD_FREQ, 1, VCFHeaderLineType.Float, GNOMAD_FREQ_DESC));
    }

    public static void addFilterHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFFilterHeaderLine(PON_GNOMAD_FILTER, "Filter Gnomad PON"));
    }
}
