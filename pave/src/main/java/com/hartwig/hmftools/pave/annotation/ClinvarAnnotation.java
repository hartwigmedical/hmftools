package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import static htsjdk.variant.vcf.VCFHeaderLineCount.UNBOUNDED;

import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.perf.StringCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class ClinvarAnnotation extends AnnotationData implements Callable<Void>
{
    private final Map<String, ClinvarChrCache> mChrCacheMap;
    private final StringCache mStringCache;
    private boolean mHasValidData;
    private final String mFilename;

    private static final String CLINVAR_VCF = "clinvar_vcf";

    public static final String CLNSIG = "CLNSIG";
    public static final String CLNSIGCONF = "CLNSIGCONF";

    public static final String CLNSIG_DESC = "Clinical significance for this single variant";
    public static final String CLNSIGCONF_DESC = "Conflicting clinical significance for this single variant";

    public ClinvarAnnotation(final ConfigBuilder configBuilder)
    {
        mChrCacheMap = Maps.newHashMap();
        mStringCache = new StringCache();

        mHasValidData = true;
        mFilename = configBuilder.getValue(CLINVAR_VCF);
    }

    @Override
    public String type() { return "Clinvar"; }

    @Override
    public boolean enabled() { return mFilename != null; }

    @Override
    public boolean hasValidData() { return mHasValidData; }

    public synchronized ClinvarChrCache getChromosomeCache(final String chromosome)
    {
        return mChrCacheMap.get(RefGenomeFunctions.stripChrPrefix(chromosome));
    }

    @Override
    public synchronized void onChromosomeComplete(final String chromosome)
    {
        ClinvarChrCache chrCache = mChrCacheMap.get(chromosome);

        if(chrCache != null)
        {
            chrCache.clear();
            mChrCacheMap.remove(chromosome);
        }
    }

    @Override
    public Void call()
    {
        if(mFilename != null)
        {
            loadEntries(mFilename);
        }

        return null;
    }

    public void loadData() { loadEntries(mFilename);}

    public static void addHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(CLNSIG, UNBOUNDED, VCFHeaderLineType.String, CLNSIG_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(CLNSIGCONF, UNBOUNDED, VCFHeaderLineType.String, CLNSIGCONF_DESC));
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(CLINVAR_VCF, false, "Clinvar annotation VCF");
    }

    private void loadEntries(final String filename)
    {
        if(filename == null)
            return;

        VcfFileReader vcfFileReader = new VcfFileReader(filename);

        if(!vcfFileReader.fileValid())
        {
            mHasValidData = false;
            return;
        }

        PV_LOGGER.debug("loading Clinvar data from file({})", filename);

        try
        {
            ClinvarChrCache currentCache = null;
            int entryCount = 0;

            for(VariantContext context : vcfFileReader.iterator())
            {
                if(context.getAlleles().size() < 2)
                    continue;

                if(!HumanChromosome.contains(context.getContig()))
                    continue;

                String chromosome = RefGenomeFunctions.stripChrPrefix(context.getContig());

                if(currentCache == null || !currentCache.Chromosome.equals(chromosome))
                {
                    currentCache = new ClinvarChrCache(chromosome, mStringCache);
                    mChrCacheMap.put(chromosome, currentCache);
                }

                int position = context.getStart();
                String ref = context.getReference().getBaseString();
                String alt = context.getAlternateAlleles().get(0).toString();

                String significance = context.getAttributeAsString(CLNSIG, "");
                String conflict = context.getAttributeAsString(CLNSIGCONF, "");

                if(significance.isEmpty() && conflict.isEmpty())
                    continue;

                currentCache.addEntry(position, ref, alt, significance, conflict);
                ++entryCount;
            }

            PV_LOGGER.info("loaded {} Clinvar entries from file({})", entryCount, filename);
        }
        catch(Exception e)
        {
            PV_LOGGER.error("failed to read Clinvar VCF file: {}", e.toString());
            mHasValidData = false;
        }
    }
}
