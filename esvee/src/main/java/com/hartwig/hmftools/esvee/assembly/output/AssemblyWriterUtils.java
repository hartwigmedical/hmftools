package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.EXTENSION;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.RemoteReadType;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;

public final class AssemblyWriterUtils
{
    protected static void addSupportHeader(final StringJoiner sj)
    {
        sj.add("SplitFrags");
        sj.add("RefSplitFrags");
        sj.add("DiscFrags");
        sj.add("RefDiscFrags");
    }

    protected static void addSupportCounts(final JunctionAssembly assembly, final StringJoiner sj)
    {
        int juncFrags = 0;
        int refSampleJuncFrags = 0;
        int discFrags = 0;
        int refSampleDiscFrags = 0;

        Set<String> uniqueReadIds = Sets.newHashSet();

        for(SupportRead read : assembly.support())
        {
            boolean isReference = read.isReference();

            if(uniqueReadIds.contains(read.id()))
                continue;

            // fragments, not reads, are being counted
            if(read.type() != EXTENSION) // let any discordant read in the fragment count
                uniqueReadIds.add(read.id());

            if(read.type().isSplitSupport())
            {
                ++juncFrags;

                if(isReference)
                    ++refSampleJuncFrags;
            }
            else
            {
                if(read.type() == DISCORDANT)
                {
                    ++discFrags;

                    if(isReference)
                        ++refSampleDiscFrags;
                }
            }
        }

        sj.add(String.valueOf(juncFrags));
        sj.add(String.valueOf(refSampleJuncFrags));
        sj.add(String.valueOf(discFrags));
        sj.add(String.valueOf(refSampleDiscFrags));
    }

    protected static void addPhasingHeader(final StringJoiner sj)
    {
        sj.add("PhaseGroupId");
        sj.add("PhaseGroupCount");

        sj.add("PhaseSetId");
        sj.add("PhaseSetCount");
        sj.add("SplitLinks");
        sj.add("FacingLinks");
        sj.add("SvType");
        sj.add("SvLength");
        sj.add("InsertedBases");
        sj.add("OverlapBases");
        sj.add("RefLinkedFrags");

        sj.add("SecondaryLinks");
    }

    private static final int NO_PHASE_ID = -1;

    protected static void addPhasingInfo(final JunctionAssembly assembly, final StringJoiner sj)
    {
        PhaseGroup phaseGroup = assembly.phaseGroup();

        int phaseGroupId = NO_PHASE_ID;
        int phaseGroupCount = 0;
        int phaseSetId = NO_PHASE_ID;
        int phaseSetCount = 0;
        String splitLinkInfo = "";
        String facingLinkInfo = "";
        String svType = "";
        int svLinkLength = 0;
        String insertedBases = "";
        String overlapBases = "";
        String linkFragmentInfo = "";
        String secondaryLinkInfo = "";

        if(phaseGroup != null)
        {
            phaseGroupId = phaseGroup.id();
            phaseGroupCount = phaseGroup.assemblyCount();

            PhaseSet phaseSet = assembly.phaseSet();

            if(phaseSet != null)
            {
                phaseSetId = phaseSet.id();
                phaseSetCount = phaseSet.assemblies().size();

                List<AssemblyLink> assemblyLinks = phaseSet.findAssemblyLinks(assembly);

                List<AssemblyLink> splitLinks = assemblyLinks.stream()
                        .filter(x -> x.type() == LinkType.SPLIT || x.type() == LinkType.INDEL).collect(Collectors.toList());

                List<AssemblyLink> facingLinks = assemblyLinks.stream().filter(x -> x.type() == LinkType.FACING).collect(Collectors.toList());
                splitLinkInfo = assemblyLinksStr(assembly, splitLinks);
                facingLinkInfo = assemblyLinksStr(assembly, facingLinks);

                if(!splitLinks.isEmpty())
                {
                    AssemblyLink svLink = splitLinks.get(0);
                    svType = svLink.svType().toString();
                    svLinkLength = svLink.length();
                    insertedBases = svLink.insertedBases();
                    overlapBases = svLink.overlapBases();
                    linkFragmentInfo = assemblyLinkFragmentStr(svLink);
                }

                List<AssemblyLink> secondarySplitLinks = phaseGroup.findSecondarySplitLinks(assembly);
                secondaryLinkInfo = assemblyLinksStr(assembly, secondarySplitLinks);
            }
        }

        sj.add(String.valueOf(phaseGroupId));
        sj.add(String.valueOf(phaseGroupCount));
        sj.add(String.valueOf(phaseSetId));
        sj.add(String.valueOf(phaseSetCount));
        sj.add(splitLinkInfo);
        sj.add(facingLinkInfo);
        sj.add(svType);
        sj.add(String.valueOf(svLinkLength));
        sj.add(insertedBases);
        sj.add(overlapBases);
        sj.add(linkFragmentInfo);
        sj.add(secondaryLinkInfo);
    }

    protected static String assemblyLinksStr(final JunctionAssembly assembly, final List<AssemblyLink> assemblyLinks)
    {
        if(assemblyLinks.isEmpty())
            return "";

        StringJoiner sj = new StringJoiner(ITEM_DELIM);

        for(AssemblyLink link : assemblyLinks)
        {
            JunctionAssembly otherAssembly = link.otherAssembly(assembly);
            sj.add(otherAssembly.junction().coords());
        }
        return sj.toString();
    }

    private static String assemblyLinkFragmentStr(final AssemblyLink assemblyLink)
    {
        int uniqueFrags = 0;
        int uniqueRefFrags = 0;
        int matchedFrags = 0;
        int matchedRefFrags = 0;

        JunctionAssembly first = assemblyLink.first();
        JunctionAssembly second = assemblyLink.second();

        Set<String> uniqueReadIds = Sets.newHashSet();

        for(SupportRead firstSupport : first.support())
        {
            if(uniqueReadIds.contains(firstSupport.id()))
                continue;

            uniqueReadIds.add(firstSupport.id());

            boolean isReference = firstSupport.isReference();

            if(second.support().stream().anyMatch(x -> x.matchesFragment(firstSupport, false)))
            {
                ++matchedFrags;

                if(isReference)
                    ++matchedRefFrags;
            }
            else
            {
                ++uniqueFrags;

                if(isReference)
                    ++uniqueRefFrags;
            }
        }

        for(SupportRead secondSupport : second.support())
        {
            if(uniqueReadIds.contains(secondSupport.id()))
                continue;

            uniqueReadIds.add(secondSupport.id());

            if(first.support().stream().noneMatch(x -> x.matchesFragment(secondSupport, true)))
            {
                ++uniqueFrags;

                if(secondSupport.isReference())
                    ++uniqueRefFrags;
            }
        }

        return format("Match=%d;MatchRef=%d;Unique=%d;UniqueRef=%d", matchedFrags, matchedRefFrags, uniqueFrags, uniqueRefFrags);
    }

    protected static String repeatsInfoStr(final List<RepeatInfo> repeats)
    {
        if(repeats.isEmpty())
            return "";

        RepeatInfo longest = null;
        RepeatInfo longestSequence = null;

        for(RepeatInfo repeat : repeats)
        {
            if(longest == null || repeat.Count > longest.Count)
                longest = repeat;

            if(longestSequence == null || repeat.length() > longest.length())
                longestSequence = repeat;
        }

        return format("%d max(%s=%d) long(%s=%d)",
                repeats.size(), longest.Bases, longest.Count, longestSequence.Bases, longestSequence.Count);
    }

    protected static void addRemoteRegionHeader(final StringJoiner sj)
    {
        sj.add("RemoteRegionCount");
        sj.add("RemoteRegionInfo");
        sj.add("RemoteRegionJuncMate");
        sj.add("RemoteRegionJuncSupps");
        sj.add("RemoteRegionDiscordant");
    }

    protected static void addRemoteRegionInfo(final JunctionAssembly assembly, final StringJoiner sj)
    {
        sj.add(String.valueOf(assembly.remoteRegions().size()));
        sj.add(remoteRegionInfoStr(assembly.remoteRegions()));

        int remoteJunctSupp = 0;
        int remoteJunctMate = 0;
        int remoteDiscordant = 0;

        for(RemoteRegion region : assembly.remoteRegions())
        {
            remoteJunctMate += region.readTypeCounts()[RemoteReadType.JUNCTION_MATE.ordinal()];
            remoteJunctSupp += region.readTypeCounts()[RemoteReadType.SUPPLEMENTARY.ordinal()];
            remoteDiscordant += region.readTypeCounts()[RemoteReadType.DISCORDANT.ordinal()];
        }

        sj.add(String.valueOf(remoteJunctMate));
        sj.add(String.valueOf(remoteJunctSupp));
        sj.add(String.valueOf(remoteDiscordant));
    }

    private static String remoteRegionInfoStr(final List<RemoteRegion> regions)
    {
        if(regions.isEmpty())
            return "";

        // log first N by read support
        Collections.sort(regions, Comparator.comparing(x -> -x.readCount()));

        StringJoiner sj = new StringJoiner(" ");

        for(int i = 0; i < min(3, regions.size()); ++i)
        {
            RemoteRegion region = regions.get(i);
            sj.add(format("%s:%d-%d=%d", region.Chromosome, region.start(), region.end(), region.readCount()));
        }

        return sj.toString();
    }

    protected static String refSideSoftClipsStr(final List<RefSideSoftClip> refSideSoftClips)
    {
        if(refSideSoftClips.isEmpty())
            return "";

        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        refSideSoftClips.forEach(x -> sj.add(format("%d:%d=%d", x.Position, x.maxLength(), x.readCount())));
        return sj.toString();
    }
}
