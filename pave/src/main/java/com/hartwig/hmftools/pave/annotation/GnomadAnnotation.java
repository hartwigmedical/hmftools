package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.StringCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.pon.GnomadChrCache;
import com.hartwig.hmftools.common.variant.pon.GnomadCommon;
import com.hartwig.hmftools.pave.VariantData;

public class GnomadAnnotation extends AnnotationData implements Callable
{
    private final Map<String,GnomadChrCache> mChrCacheMap;
    private final RefGenomeVersion mRefGenomeVersion;
    private final Map<String,String> mChromosomeFiles;
    private boolean mHasValidData;
    private final boolean mEnabled;
    private final String mGnomadFilename;
    private final StringCache mStringCache;
    private final boolean mNoFilter;

    public static final String GNOMAD_FREQUENCY_FILE = "gnomad_freq_file";
    public static final String GNOMAD_FREQUENCY_DIR = "gnomad_freq_dir";
    public static final String GNOMAD_NO_FILTER = "gnomad_no_filter";

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

        mNoFilter = configBuilder.hasFlag(GNOMAD_NO_FILTER);
    }

    @Override
    public String type() { return "Gnomad frequency"; }

    @Override
    public boolean enabled() { return mGnomadFilename != null || !mChromosomeFiles.isEmpty(); }

    public boolean applyFilter() { return !mNoFilter; }

    @Override
    public boolean hasValidData() { return mHasValidData; }

    public void annotateVariant(final VariantData variant, final GnomadChrCache chrCache)
    {
        Double gnomadFreq = chrCache.getFrequency(variant.isMnv(), variant.Ref, variant.Alt, variant.Position);

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

        mHasValidData = GnomadCommon.loadChromosomeEntries(filename, fileChromosome, mStringCache, mChrCacheMap, mChromosomeFiles);
    }

    private void loadAllFrequencyFiles(final String gnomadDir)
    {
        mHasValidData = GnomadCommon.loadAllFrequencyFiles(gnomadDir, mRefGenomeVersion, mChromosomeFiles);
    }
}
