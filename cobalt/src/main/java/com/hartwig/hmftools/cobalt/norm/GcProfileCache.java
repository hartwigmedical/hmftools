package com.hartwig.hmftools.cobalt.norm;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.loadChrGcProfileMap;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.gc.GCProfile;

public class GcProfileCache
{
    private final Map<String, List<GCProfile>> mChrGcProfiles;

    private String mCurentChromosome;
    private int mPositionIndex;

    public GcProfileCache(final String filename)
    {
        mChrGcProfiles = Maps.newHashMap();
        mPositionIndex = 0;
        mCurentChromosome = "";

        if(filename == null)
        {
            return;
        }

        try
        {
            CB_LOGGER.info("loading GC profile: {}", filename);
            mChrGcProfiles.putAll(loadChrGcProfileMap(filename));
        }
        catch(IOException e)
        {
            CB_LOGGER.error("error loading GC profile: {}", e.toString());
        }
    }

    public Map<String, List<GCProfile>> chrGcProfiles()
    {
        return mChrGcProfiles;
    }

    public GCProfile findGcProfile(final String chromosome, final int position)
    {
        if(!mCurentChromosome.equals(chromosome))
        {
            mPositionIndex = 0;
            mCurentChromosome = chromosome;
        }

        List<GCProfile> profiles = mChrGcProfiles.get(chromosome);

        int index = mPositionIndex;
        while(index < profiles.size())
        {
            GCProfile profile = profiles.get(index);

            if(profile.start() == position)
            {
                mPositionIndex = index;
                return profile;
            }
            else if(position < profile.start())
            {
                return null;
            }

            ++index;
        }

        return null;
    }
}
