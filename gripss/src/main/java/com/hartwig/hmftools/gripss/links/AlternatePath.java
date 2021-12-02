package com.hartwig.hmftools.gripss.links;

import static com.hartwig.hmftools.gripss.links.Link.LINK_TYPE_PAIR;
import static com.hartwig.hmftools.gripss.links.TransitiveLink.TRANS_LINK_PREFIX;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.gripss.common.Breakend;

public class AlternatePath
{
    // an alternative path is modelled as the pair of breakends from the same SV, in a specific order
    // and the set of links made by other SVs that start and end at the same locations
    public final Breakend First;
    public final Breakend Second;
    public final List<Link> Links;

    public AlternatePath(final Breakend first, final Breakend second, final List<Link> links)
    {
        First = first;
        Second = second;
        Links = links;
    }

    public List<Link> transitiveLinks()
    {
        return Links.stream().filter(x -> x.Id.startsWith(TRANS_LINK_PREFIX)).collect(Collectors.toList());
    }

    public String pathString()
    {
        StringJoiner sj = new StringJoiner("");

        for(int i = 0; i < Links.size(); ++i)
        {
            Link link = Links.get(i);

            if (i == 0)
                sj.add(link.breakendStart().VcfId);

            if (link.Id.equals(LINK_TYPE_PAIR))
                sj.add("-");
            else
                sj.add(String.format("<%s>", link.Id));

            sj.add(link.breakendEnd().VcfId);
        }

        return sj.toString();
    }

    public String toString()
    {
        return String.format("breaks(%s - %s) links(%d)", First, Second, Links.size());
    }

}
