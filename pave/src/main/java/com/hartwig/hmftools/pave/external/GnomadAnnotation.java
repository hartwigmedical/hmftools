package com.hartwig.hmftools.pave.external;

import static java.lang.Math.max;

import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.pave.VariantData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class GnomadAnnotation
{
    private final Map<String,Map<Integer,List<GnomadVariant>>> mFrequencies;

    public static final String GNOMAD_FREQUENCY_FILE = "gnomad_freq_file";
    public static final String GNOMAD_VCF_TAG = "GND_FREQ";

    public GnomadAnnotation(final CommandLine cmd)
    {
        mFrequencies = Maps.newHashMap();

        if(cmd != null)
        {
            loadFrequency(cmd.getOptionValue(GNOMAD_FREQUENCY_FILE));
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

    private void loadFrequency(final String filename)
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

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(",", -1);

                String chromosome = values[0];
                int position = Integer.parseInt(values[1]);
                String ref = values[2];
                String alt = values[3];
                double frequency = Double.parseDouble(values[4]);

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

            PV_LOGGER.info("loaded {} gnoamd frequency records from file", itemCount, filename);
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load gnoamd frequency file({}): {}", filename, e.toString());
            return;
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(GNOMAD_FREQUENCY_FILE, true, "Optional: write all batch-run output files");
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

    /*
    private class GnomadVariantFrequencies
    {
        public final char[] RefAltKey;
        public final double Frequency;

        public GnomadVariant(final String ref, final String alt, final double frequency)
        {
            int keyLength = ref.length() + alt.length() + 2;
            RefAltKey = new char[keyLength];

            int index = 0;
            for(int i = 0; i < ref.length(); ++i)
            {
                RefAltKey[index++] = ref.charAt(i);
            }

            RefAltKey[index++] = '>';

            for(int i = 0; i < alt.length(); ++i)
            {
                RefAltKey[index++] = alt.charAt(i);
            }

            Frequency = frequency;
        }

        public
    }
    */
}
