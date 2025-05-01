package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.List;
import java.util.NoSuchElementException;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.TaskQueue;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;

public class PonThread extends Thread
{
    private final PonConfig mConfig;
    private final TaskQueue mRegions;
    private final PonWriter mPonWriter;

    private final List<String> mSampleVcfs;
    private final ClinvarAnnotation mClinvarAnnotation;
    private final HotspotCache mHotspotCache;
    private final EnsemblDataCache mEnsemblDataCache;

    public PonThread(
            final PonConfig config, final List<String> sampleVcfs, final TaskQueue taskQueue, final PonWriter ponWriter,
            final ClinvarAnnotation clinvarAnnotation, final HotspotCache hotspotCache, final EnsemblDataCache ensemblDataCache)
    {
        mRegions = taskQueue;
        mConfig = config;
        mPonWriter = ponWriter;
        mSampleVcfs = sampleVcfs;
        mEnsemblDataCache = ensemblDataCache;
        mHotspotCache = hotspotCache;
        mClinvarAnnotation = clinvarAnnotation;
    }

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

                mPonWriter.onVariantsComplete(region, regionPonTask.variants());
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
