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

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.StringCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class GnomadCommon
{
    public static final String GNOMAD_FREQUENCY_FILE = "gnomad_freq_file";
    public static final String GNOMAD_FREQUENCY_DIR = "gnomad_freq_dir";
    public static final String GNOMAD_NO_FILTER = "gnomad_no_filter";

    public static final String PON_GNOMAD_FILTER = "PONGnomad";

    public static final String GNOMAD_FILE_ID = "gnomad_variants";

    protected static final Logger LOGGER = LogManager.getLogger(GnomadCommon.class);

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

    public static boolean loadChromosomeEntries(
            final String filename, final String fileChromosome, final StringCache stringCache, final Map<String,GnomadChrCache> chrCacheMap,
            final Map<String,String> chromosomeFiles)
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
                    currentChrCache = new GnomadChrCache(chromosome, stringCache);
                    chrCacheMap.put(chromosome, currentChrCache);
                }

                currentChrCache.addEntry(position, ref, alt,frequency);

                ++itemCount;
            }

            if(chromosomeFiles.isEmpty())
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

    public static boolean loadAllFrequencyFiles(
            final String gnomadDir, final RefGenomeVersion refGenomeVersion, final Map<String,String> chromosomeFiles)
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

                String chrStr = refGenomeVersion.versionedChromosome(humanChr.toString());

                if(chrFile == null)
                {
                    LOGGER.error("missing Gnomad chromosome({}) file", chrStr);
                    return false;
                }

                chromosomeFiles.put(chrStr, chrFile);
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
        configBuilder.addFlag(GNOMAD_NO_FILTER, "No Gnomad filter is applied");
    }

    public static void addHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(GNOMAD_FREQ, 1, VCFHeaderLineType.Float, GNOMAD_FREQ_DESC));
        header.addMetaDataLine(new VCFFilterHeaderLine(PON_GNOMAD_FILTER, "Filter Gnomad PON"));
    }
}
