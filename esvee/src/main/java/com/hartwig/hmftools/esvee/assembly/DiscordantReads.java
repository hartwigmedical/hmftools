package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConstants.READ_SOFT_CLIP_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.common.DiscordantGroup;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

public final class DiscordantReads
{
    public static List<DiscordantGroup> buildFromDiscordantReads(final Junction junction, final List<Read> rawReads)
    {
        List<DiscordantGroup> groups = Lists.newArrayList();

        for(Read read : rawReads)
        {
            if(!isDiscordant(read))
                continue;

            if(read.orientation() != junction.Orientation)
                continue;

            if(junction.isForward())
            {
                if(read.alignmentEnd() > junction.position() + READ_SOFT_CLIP_JUNCTION_BUFFER)
                    continue;
            }
            else
            {
                if(read.alignmentStart() < junction.position() - READ_SOFT_CLIP_JUNCTION_BUFFER)
                    continue;
            }

            DiscordantGroup matchedGroup = groups.stream().filter(x -> x.matches(read)).findFirst().orElse(null);

            if(matchedGroup != null)
            {
                matchedGroup.addRead(read);
            }
            else
            {
                groups.add(new DiscordantGroup(junction, read));
            }
        }

        return groups;
    }
}
