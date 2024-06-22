package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ_DESC;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.resources.GnomadCacheBuilder.GNOMAD_FILE_ID;
import static com.hartwig.hmftools.pave.resources.GnomadCacheBuilder.formFileId;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.StringCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.pave.VariantData;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GnomadAnnotation extends AnnotationData implements Callable
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
    private static final String GNOMAD_PON_FILTER = "gnomad_pon_filter";

    public static final String PON_GNOMAD_FILTER = "PONGnomad";

    public GnomadAnnotation(final ConfigBuilder configBuilder)
    {
        mChrCacheMap = Maps.newHashMap();
        mChromosomeFiles = Maps.newHashMap();
        mHasValidData = true;
        mStringCache = new StringCache();

        mRefGenomeVersion = RefGenomeVersion.from(configBuilder);

        if(configBuilder.hasValue(GNOMAD_FREQUENCY_FILE))
        {
            mEnabled = true;
            mGnomadFilename = configBuilder.getValue(GNOMAD_FREQUENCY_FILE);
        }
        else if(configBuilder.hasValue(GNOMAD_FREQUENCY_DIR))
        {
            mEnabled = true;
            mGnomadFilename = null;
            String gnomadDir = configBuilder.getValue(GNOMAD_FREQUENCY_DIR);
            loadAllFrequencyFiles(gnomadDir);
        }
        else
        {
            mGnomadFilename = null;
            mEnabled = false;
        }
    }

    @Override
    public String type() { return "Gnomad frequency"; }

    @Override
    public boolean enabled() { return mGnomadFilename != null || !mChromosomeFiles.isEmpty(); }

    @Override
    public boolean hasValidData() { return mHasValidData; }

    public void annotateVariant(final VariantData variant, final GnomadChrCache chrCache)
    {
        Double gnomadFreq = chrCache.getFrequency(variant);

        if(gnomadFreq != null)
        {
            variant.setGnomadFrequency(gnomadFreq);
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
            PV_LOGGER.error("missing Gnomad file for chromosome({})", chromosome);
            System.exit(1);
        }

        loadChromosomeEntries(chrFilename, chromosome);
        return mChrCacheMap.get(chromosome);
    }

    @Override
    public synchronized void onChromosomeComplete(final String chromosome)
    {
        GnomadChrCache chrCache = mChrCacheMap.get(chromosome);

        if(chrCache != null)
        {
            chrCache.clear();
            mChrCacheMap.remove(chromosome);
        }
    }

    @Override
    public Long call()
    {
        if(mGnomadFilename != null)
        {
            loadChromosomeEntries(mGnomadFilename, null);
        }
        else if(!mChromosomeFiles.isEmpty())
        {
            for(String chromosome : mInitialChromosomes)
            {
                getChromosomeCache(chromosome);
            }
        }

        return (long)0;
    }

    private void loadChromosomeEntries(final String filename, final String fileChromosome)
    {
        // if file chromosome is supplied then it is not read from the input file
        if(filename == null)
            return;

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
                PV_LOGGER.info("loaded {} Gnomad frequency records from file({})", itemCount, filename);
            }
            else if(currentChrCache != null)
            {
                PV_LOGGER.debug("chr({}) loaded {} Gnomad frequency records",
                        currentChrCache.Chromosome, currentChrCache.entryCount());
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load Gnomad frequency file({}): {}", filename, e.toString());
            mHasValidData = false;
        }
    }

    private void loadAllFrequencyFiles(final String gnomadDir)
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
                    PV_LOGGER.error("missing Gnomad chromosome({}) file", chrStr);
                    mHasValidData = false;
                    return;
                }

                mChromosomeFiles.put(chrStr, chrFile);
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to find Gnomad chromosome files in dir({}): {}", gnomadDir, e.toString());
        }
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(GNOMAD_FREQUENCY_FILE, false, "Gnomad frequency file");
        configBuilder.addPath(GNOMAD_FREQUENCY_DIR, false, "Gnomad frequency directory");
    }

    public static void addHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(GNOMAD_FREQ, 1, VCFHeaderLineType.Float, GNOMAD_FREQ_DESC));
        header.addMetaDataLine(new VCFFilterHeaderLine(PON_GNOMAD_FILTER, "Filter Gnomad PON"));
    }
}
