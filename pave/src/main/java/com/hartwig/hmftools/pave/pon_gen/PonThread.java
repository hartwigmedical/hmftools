package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.List;
import java.util.NoSuchElementException;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.TaskQueue;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;

public class PonThread extends Thread
{
    private final PonConfig mConfig;
    private final TaskQueue mRegions;

    private final List<String> mSampleVcfs;
    private final ClinvarAnnotation mClinvarAnnotation;
    private final HotspotCache mHotspotCache;
    private final EnsemblDataCache mEnsemblDataCache;

    private final List<VariantPonData> mFilteredVariants;

    public PonThread(
            final PonConfig config, final List<String> sampleVcfs, final TaskQueue taskQueue,
            final ClinvarAnnotation clinvarAnnotation, final HotspotCache hotspotCache, final EnsemblDataCache ensemblDataCache)
    {
        mRegions = taskQueue;
        mConfig = config;
        mSampleVcfs = sampleVcfs;
        mEnsemblDataCache = ensemblDataCache;
        mHotspotCache = hotspotCache;
        mClinvarAnnotation = clinvarAnnotation;

        mFilteredVariants = Lists.newArrayList();
    }

    public List<VariantPonData> filteredVariants() { return mFilteredVariants; }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                ChrBaseRegion region = (ChrBaseRegion)mRegions.removeItem();

                RegionPonTask regionPonTask = new RegionPonTask(
                        mConfig, region, mSampleVcfs, mClinvarAnnotation, mHotspotCache, mEnsemblDataCache);

                regionPonTask.run();

                mFilteredVariants.addAll(regionPonTask.variants());
            }
            catch(NoSuchElementException e)
            {
                PV_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }
}
