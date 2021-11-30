package com.hartwig.hmftools.gripss.links;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.gripss.common.Breakend;

public class AlternatePath
{
    public final Breakend First;
    public final Breakend Second;
    public final List<Link> Links;

    public AlternatePath(final Breakend first, final Breakend second, final List<Link> links)
    {
        First = first;
        Second = second;
        Links = links;
    }

    public List<String> pathVcfIds()
    {
        List<String> pathStrings = Links.stream().map(x -> x.breakendStart().VcfId).collect(Collectors.toList());
        pathStrings.add(Links.get(Links.size() - 1).breakendEnd().VcfId);
        return pathStrings;
    }

    public List<Link> transitiveLinks()
    {
        // TODO: either make 'trs' a constant or add a link type
        return Links.stream().filter(x -> x.Id.startsWith("trs")).collect(Collectors.toList());
    }

    public String pathString()
    {
        StringJoiner sj = new StringJoiner("");

        for(int i = 0; i < Links.size(); ++i)
        {
            Link link = Links.get(i);

            if (i == 0)
                sj.add(link.breakendStart().VcfId);

            if (link.Id.equals("PAIR"))
                sj.add("-");
            else
                sj.add(String.format("<%s>", link.Id));

            sj.add(link.breakendEnd().VcfId);
        }

        return sj.toString();
    }

    public static LinkStore createLinkStore(final List<AlternatePath> alternatePaths)
    {
        LinkStore linkStore = new LinkStore();

        for(AlternatePath altPath : alternatePaths)
        {
            List<Link> transLinks = altPath.transitiveLinks();
            transLinks.forEach(x -> linkStore.addLink(altPath.First, x));
        }

        return linkStore;
    }

}
