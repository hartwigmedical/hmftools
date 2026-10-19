package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.VariantHotspotFile.SOMATIC_HOTSPOT_FILE;
import static com.hartwig.hmftools.common.variant.VariantHotspotFile.SOMATIC_HOTSPOT_FILE_DESC;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_FILE;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_FILE_DESC;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.APP_NAME;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VariantHotspotFile;

import org.jetbrains.annotations.NotNull;

public class PonHotspotChecker
{
    private static final String OUTPUT_PON_FILE = "output_pon_file";

    public PonHotspotChecker(ConfigBuilder configBuilder)
    {
        String ponFile = configBuilder.getValue(PON_FILE);
        String hotspotFile = configBuilder.getValue(SOMATIC_HOTSPOT_FILE);
        String outputPonFile = configBuilder.getValue(OUTPUT_PON_FILE);

        PV_LOGGER.info("checking PON file({}) for hotspot matches from {}", ponFile, hotspotFile);

        List<String> validPonLines = Lists.newArrayList();
        boolean hasHotspotMatches = false;

        try
        {
            ListMultimap<Chromosome, SimpleVariant> hotspotMap = VariantHotspotFile.loadHotspotVcf(hotspotFile);

            BufferedReader fileReader = createBufferedReader(ponFile);

            String line = fileReader.readLine();
            String header = line;

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, TSV_DELIM);
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int positionIndex = fieldsIndexMap.get(FLD_POSITION);
            int refIndex = fieldsIndexMap.get(FLD_REF);
            int altIndex = fieldsIndexMap.get(FLD_ALT);

            List<SimpleVariant> chrHotspots = null;
            String currentChromosome = "";

            // Chromosome      Position        Ref     Alt     SampleCount     MaxSampleReadCount      TotalReadCount  etc
            while((line = fileReader.readLine()) != null)
            {
                String[] values = line.split(TSV_DELIM, -1);

                String chromosome = values[chrIndex];

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    chrHotspots = hotspotMap.get(HumanChromosome.fromString(chromosome));
                }

                int position = Integer.parseInt(values[positionIndex]);
                String ref = values[refIndex];
                String alt = values[altIndex];

                // check for a hotspot match
                if(chrHotspots.stream().anyMatch(x -> x.matches(chromosome, position, ref, alt)))
                {
                    PV_LOGGER.debug("PON entry({}:{}:{}:{}) matches hotspot", chromosome, position, ref, alt);
                    hasHotspotMatches = true;
                    continue;
                }

                if(outputPonFile != null)
                    validPonLines.add(line);
            }

            if(hasHotspotMatches && outputPonFile != null)
            {
                PV_LOGGER.info("writing filtered PON file: {}", outputPonFile);

                BufferedWriter ponWriter = createBufferedWriter(outputPonFile);

                ponWriter.write(header);
                ponWriter.newLine();

                for(String ponLine : validPonLines)
                {
                    ponWriter.write(ponLine);
                    ponWriter.newLine();
                }

                ponWriter.close();
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load hotspot file: {}", e.toString());
            System.exit(1);
        }

        PV_LOGGER.info("PON hotspot check complete");
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PonHotspotChecker ponHotspotChecker = new PonHotspotChecker(configBuilder);
    }

    private static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(OUTPUT_PON_FILE, false, "Output PON with hotspots filtered out");
        configBuilder.addPath(PON_FILE, true, PON_FILE_DESC);
        configBuilder.addPath(SOMATIC_HOTSPOT_FILE, false, SOMATIC_HOTSPOT_FILE_DESC);
        addLoggingOptions(configBuilder);
    }
}
