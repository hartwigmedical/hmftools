package com.hartwig.hmftools.common.variant.pon;

import static java.lang.String.format;

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

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.StringCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class GnomadChrCache
{
    public final String Chromosome;

    private final Map<Integer,List<GnomadVariant>> mFrequencies;
    private final StringCache mStringCache;

    public static final String GNOMAD_FREQUENCY_FILE = "gnomad_freq_file";
    public static final String GNOMAD_FREQUENCY_DIR = "gnomad_freq_dir";
    public static final String GNOMAD_NO_FILTER = "gnomad_no_filter";

    public static final String PON_GNOMAD_FILTER = "PONGnomad";

    protected static final Logger LOGGER = LogManager.getLogger(GnomadCommon.class);


    public GnomadChrCache(final String chromosome, final StringCache stringCache)
    {
        Chromosome = chromosome;
        mFrequencies = Maps.newHashMap();
        mStringCache = stringCache;
    }

    public void addEntry(final int position, final String ref, final String alt, final double frequency)
    {
        List<GnomadVariant> posList = mFrequencies.get(position);

        if(posList == null)
        {
            posList = Lists.newArrayList();
            mFrequencies.put(position, posList);
        }

        posList.add(new GnomadVariant(mStringCache.intern(ref), mStringCache.intern(alt), frequency));
    }

    public void clear() { mFrequencies.clear(); }
    public int entryCount() { return mFrequencies.size(); }

    public String toString() { return format("chr(%s) entries(%d)", Chromosome, mFrequencies.size()); }

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

    public Double getFrequency(final boolean isMnv, final String ref, final String alt, int position)
    {
        if(isMnv)
        {
            // check for successive bases for an MNV
            double minFreq = 0;
            for(int i = 0; i < ref.length(); ++i)
            {
                Double freq = getFrequency(position + i, ref.substring(i, i + 1), alt.substring(i, i + 1));

                if(freq == null)
                    return null;

                if(minFreq == 0 ||freq < minFreq)
                    minFreq = freq;
            }

            return minFreq;
        }
        else
        {
            return getFrequency(position, ref, alt);
        }
    }

    public Double getFrequency(int position, final String ref, final String alt)
    {
        List<GnomadVariant> posList = mFrequencies.get(position);

        if(posList == null)
            return null;

        GnomadVariant match = posList.stream().filter(x -> x.matches(ref, alt)).findFirst().orElse(null);
        return match != null ? match.Frequency : null;
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
