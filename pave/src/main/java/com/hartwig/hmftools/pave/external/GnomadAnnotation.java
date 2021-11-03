package com.hartwig.hmftools.pave.external;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.external.GnomadParser.GNOMAD_FILE_ID;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.pave.VariantData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class GnomadAnnotation
{
    private final Map<String,Map<Integer,List<GnomadVariant>>> mFrequencies;

    public static final String GNOMAD_FREQUENCY_FILE = "gnomad_freq_file";
    public static final String GNOMAD_FREQUENCY_DIR = "gnomad_freq_dir";
    public static final String GNOMAD_VCF_TAG = "GND_FREQ";

    public GnomadAnnotation(final CommandLine cmd)
    {
        mFrequencies = Maps.newHashMap();

        if(cmd != null)
        {
            if(cmd.hasOption(GNOMAD_FREQUENCY_FILE))
            {
                loadFrequency(cmd.getOptionValue(GNOMAD_FREQUENCY_FILE), null);
            }
            else if(cmd.hasOption(GNOMAD_FREQUENCY_DIR))
            {
                String gnomadDir = cmd.getOptionValue(GNOMAD_FREQUENCY_DIR);
                RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));

                try
                {
                    final Stream<Path> stream = Files.walk(Paths.get(gnomadDir), 1, FileVisitOption.FOLLOW_LINKS);

                    List<String> files = stream.filter(x -> !x.toFile().isDirectory())
                            .map(x -> x.toFile().toString())
                            .filter(x -> x.contains(GNOMAD_FILE_ID))
                            .collect(Collectors.toList());

                    for(HumanChromosome humanChr : HumanChromosome.values())
                    {
                        String chrStr = refGenomeVersion.versionedChromosome(humanChr.toString());
                        String fileChrStr = "chr" + humanChr.toString();

                        String chrFile = files.stream().filter(x -> x.contains(fileChrStr)).findFirst().orElse(null);

                        if(chrFile == null)
                        {
                            PV_LOGGER.error("missing Gnomad chromosome({}) file", chrStr);
                            continue;
                        }

                        loadFrequency(chrFile, chrStr);
                    }
                }
                catch(IOException e)
                {
                    PV_LOGGER.error("failed to find gnoamd chromosome files in dir({}): {}", gnomadDir, e.toString());
                }


            }
        }
    }

    public boolean hasData() { return !mFrequencies.isEmpty(); }

    public Double getFrequency(final VariantData variant)
    {
        Map<Integer,List<GnomadVariant>> posMap = mFrequencies.get(variant.Chromosome);

        if(posMap == null)
            return null;

        List<GnomadVariant> posList = posMap.get(variant.Position);

        if(posList == null)
            return null;

        GnomadVariant match = posList.stream().filter(x -> x.matches(variant.Ref, variant.Alt)).findFirst().orElse(null);
        return match != null ? match.Frequency : null;
    }

    private void loadFrequency(final String filename, final String fileChromosome)
    {
        if(filename == null)
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            int itemCount = 0;
            String line = fileReader.readLine(); // skip header
            String currentChr = "";
            int currentPos = 0;
            int maxAtPos = 0;

            Map<Integer,List<GnomadVariant>> posMap = null;
            List<GnomadVariant> posList = null;

            Integer chrIndex = fileChromosome != null ? -1 : 0;
            int index = fileChromosome != null ? 0 : 1;
            int posIndex = index++;
            int refIndex = index++;
            int altIndex = index++;
            int freqIndex = index++;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(",", -1);

                String chromosome = fileChromosome != null ? fileChromosome : values[chrIndex];
                int position = Integer.parseInt(values[posIndex]);
                String ref = values[refIndex];
                String alt = values[altIndex];
                double frequency = Double.parseDouble(values[freqIndex]);

                if(!chromosome.equals(currentChr))
                {
                    currentChr = chromosome;
                    currentPos = position;
                    posMap = Maps.newHashMap();
                    mFrequencies.put(chromosome, posMap);

                    posList = Lists.newArrayList();
                    posMap.put(position, posList);
                }
                else if(currentPos != position)
                {
                    currentPos = position;
                    posList = Lists.newArrayList();
                    posMap.put(position, posList);
                }

                posList.add(new GnomadVariant(ref, alt, frequency));
                maxAtPos = max(maxAtPos, posList.size());

                ++itemCount;
            }

            PV_LOGGER.info("loaded {} gnomad frequency records from file({})", itemCount, filename);
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load gnoamd frequency file({}): {}", filename, e.toString());
            return;
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(GNOMAD_FREQUENCY_FILE, true, "Gnomad frequency file");
        options.addOption(GNOMAD_FREQUENCY_DIR, true, "Gnomad frequency directory");
    }

    private class GnomadVariant
    {
        public final String Ref;
        public final String Alt;
        public final double Frequency;

        public GnomadVariant(final String ref, final String alt, final double frequency)
        {
            Ref = ref;
            Alt = alt;
            Frequency = frequency;
        }

        public boolean matches(final String ref, final String alt) { return ref.equals(Ref) && alt.equals(Alt); }
    }
}
