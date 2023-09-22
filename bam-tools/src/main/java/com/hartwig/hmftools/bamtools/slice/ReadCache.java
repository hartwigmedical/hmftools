package com.hartwig.hmftools.bamtools.slice;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class ReadCache
{
    private final SliceConfig mConfig;
    private final Map<String,List<RemotePosition>> mChrRemotePositions;

    public ReadCache(final SliceConfig config)
    {
        mConfig = config;
        mChrRemotePositions = Maps.newHashMap();
    }

    public Map<String,List<RemotePosition>> chrRemotePositions() { return mChrRemotePositions; }

    public synchronized void addRemotePosition(final RemotePosition remotePosition)
    {
        // positions already part of the initial slice are assumed to have been checked and omitted
        List<RemotePosition> positions = mChrRemotePositions.get(remotePosition.Chromosome);

        if(positions == null)
        {
            mChrRemotePositions.put(remotePosition.Chromosome, Lists.newArrayList(remotePosition));
        }
        else
        {
            int index = 0;
            while(index < positions.size())
            {
                if(remotePosition.Position < positions.get(index).Position)
                    break;

                if(remotePosition.Position == positions.get(index).Position
                && remotePosition.ReadId.equals(positions.get(index).ReadId))
                {
                    return;
                }

                ++index;
            }

            positions.add(index, remotePosition);
        }
    }
}
