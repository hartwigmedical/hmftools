package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isDiscordantFragment;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.types.DiscordantGroup;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class DiscordantReads
{
    private List<DiscordantGroup> mGroups;

    public DiscordantReads()
    {
        mGroups = Lists.newArrayList();
    }

    public List<DiscordantGroup> groups() { return mGroups; }

    public void processReads(final Junction junction, final List<Read> rawReads)
    {
        DiscordantGroup currentGroup = null;

        for(Read read : rawReads)
        {
            if(read.isUnmapped() || !read.isPairedRead() || read.isMateUnmapped() || !isDiscordantFragment(read))
                continue;

            // search amongst existing groups and match on orientation and remote region
            if(currentGroup != null && currentGroup.matches(read))
            {
                currentGroup.addRead(read);
                continue;
            }

            currentGroup = mGroups.stream().filter(x -> x.matches(read)).findFirst().orElse(null);

            if(currentGroup != null)
            {
                currentGroup.addRead(read);
            }
            else
            {
                currentGroup = new DiscordantGroup(junction, read);
                mGroups.add(currentGroup);
            }
        }
    }

    public void mergeGroups()
    {
        int index = 0;
        while(index < mGroups.size() - 1)
        {
            DiscordantGroup first = mGroups.get(index);
            boolean foundMerge = false;

            int nextIndex = index + 1;
            while(nextIndex < mGroups.size())
            {
                DiscordantGroup next = mGroups.get(nextIndex);

                if(first.matches(next))
                {
                    mGroups.remove(nextIndex);
                    next.reads().stream().filter( x -> !first.hasRead(x)).forEach(x -> first.addRead(x));
                    foundMerge = true;
                }
                else
                {
                    ++nextIndex;
                }
            }

            if(!foundMerge)
                ++index;
        }
    }
}
