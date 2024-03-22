package com.hartwig.hmftools.esvee.output;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.convertedIndelCrossesJunction;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_DISCORDANT_READ;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_JUNCTION_SUPP;
import static com.hartwig.hmftools.esvee.common.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.common.SupportType.INDEL;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.common.AssemblyLink;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.LinkType;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.PhaseSet;
import com.hartwig.hmftools.esvee.common.RefSideSoftClip;
import com.hartwig.hmftools.esvee.common.RemoteRegion;
import com.hartwig.hmftools.esvee.common.RepeatInfo;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadUtils;

public final class AssemblyWriterUtils
{
    protected static void addSupportHeader(final StringJoiner sj)
    {
        sj.add("SplitFrags");
        sj.add("RefSplitFrags");
        sj.add("DiscFrags");
        sj.add("RefDiscFrags");

        sj.add("JuncSupps");
        sj.add("JuncIndels");
        sj.add("JuncMateConcord");
        sj.add("JuncMateRefSide");
        sj.add("JuncMateDiscordRemote");
        sj.add("JuncMateDiscordRefSide");
        sj.add("JuncMateUnmappedRemote");
        sj.add("JuncMateUnmappedRefSide");
    }

    protected static void addSupportCounts(final JunctionAssembly assembly, final StringJoiner sj)
    {
        int juncFrags = 0;
        int refSampleJuncFrags = 0;
        int discFrags = 0;
        int refSampleDiscFrags = 0;

        int juncSupps = 0;
        int juncMateConcordant = 0;
        int juncMateDiscordantRemote = 0;
        int juncMateDiscordantRefSide = 0;
        int juncMateRefSide = 0;
        int juncMateUnmappedRemote = 0;
        int juncMateUnmappedRefSide = 0;
        int indelReads = 0;

        int juncUnlinkedMates = 0;
        int juncUnlinkedSupps = 0;
        int discUnlinkedMates = 0;

        Set<String> uniqueReadIds = Sets.newHashSet();

        for(AssemblySupport support : assembly.support())
        {
            boolean isReference = support.read().isReference();
            Read read = support.read();

            if(support.type().isSplitSupport())
            {
                if(!uniqueReadIds.contains(read.getName()))
                {
                    uniqueReadIds.add(read.getName()); // since fragments are being counted

                    ++juncFrags;

                    if(isReference)
                        ++refSampleJuncFrags;
                }

                if(support.type() == INDEL || convertedIndelCrossesJunction(assembly, read))
                    ++indelReads;

                if(read.isSupplementary())
                {
                    ++juncSupps;

                    if(read.supplementaryRead() == null)
                        ++juncUnlinkedSupps;
                }
                else
                {
                    // don't attempt to set mate info for supplementaries - these counts are for the primary reads (ie the fragment itself)
                    boolean matePastJunction = (read.orientation() == POS_ORIENT) == assembly.isForwardJunction();

                    if(read.isMateUnmapped())
                    {
                        if(matePastJunction)
                            ++juncMateUnmappedRemote;
                        else
                            ++juncMateUnmappedRefSide;
                    }
                    else
                    {
                        if(isDiscordant(read))
                        {
                            if(matePastJunction)
                                ++juncMateDiscordantRemote;
                            else
                                ++juncMateDiscordantRefSide;
                        }
                        else
                        {
                            // check if the mate

                            if((assembly.isForwardJunction() && read.mateAlignmentStart() > assembly.junction().Position)
                                    || (!assembly.isForwardJunction() && read.mateAlignmentEnd() < assembly.junction().Position))
                            {
                                ++juncMateConcordant;
                            }
                            else
                            {
                                ++juncMateRefSide;
                            }
                        }
                    }

                    if(read.isMateMapped() && !read.hasMateSet())
                        ++juncUnlinkedMates;
                }
            }
            else
            {
                if(support.type() == DISCORDANT)
                {
                    ++discFrags;

                    if(isReference)
                        ++refSampleDiscFrags;

                    if(!read.isSupplementary() && !read.hasMateSet())
                        ++discUnlinkedMates;
                }
            }
        }

        sj.add(String.valueOf(juncFrags));
        sj.add(String.valueOf(refSampleJuncFrags));
        sj.add(String.valueOf(discFrags));
        sj.add(String.valueOf(refSampleDiscFrags));

        sj.add(String.valueOf(juncSupps));
        sj.add(String.valueOf(indelReads));

        sj.add(String.valueOf(juncMateConcordant));
        sj.add(String.valueOf(juncMateRefSide));
        sj.add(String.valueOf(juncMateDiscordantRemote));
        sj.add(String.valueOf(juncMateDiscordantRefSide));
        sj.add(String.valueOf(juncMateUnmappedRemote));
        sj.add(String.valueOf(juncMateUnmappedRefSide));

        /*
        sj.add(String.valueOf(juncUnlinkedMates));
        sj.add(String.valueOf(juncUnlinkedSupps));
        sj.add(String.valueOf(discUnlinkedMates));
            */

    }

    protected static void addReadStatsHeader(final StringJoiner sj)
    {
        sj.add("SoftClipMatches");
        sj.add("SoftClipMismatches");
        sj.add("SoftClip2ndMaxLength");
        sj.add("RefBaseMismatches");
        sj.add("BaseTrimCount");

        sj.add("AvgNmCount");
        sj.add("AvgIndelLength");
        sj.add("AvgBaseQual");
        sj.add("AvgMapQual");
        sj.add("MismatchReads");
    }

    protected static void addReadStats(final JunctionAssembly assembly, final StringJoiner sj)
    {
        int nmCountTotal = 0;
        int indelLengthTotal = 0;
        int baseQualTotal = 0;
        int mapQualTotal = 0;
        int baseTrimTotal = 0;
        int softClipMatchTotal = 0;
        int softClipMismatchTotal = 0;
        int refBaseMismatchTotal = 0;
        int softClipSecondMaxLength = 0;

        int maxExtBaseMatchCount = 0;

        for(AssemblySupport support : assembly.support())
        {
            Read read = support.read();
            nmCountTotal += read.numberOfEvents();
            mapQualTotal += read.mappingQuality();
            baseQualTotal += ReadUtils.avgBaseQuality(read);
            indelLengthTotal += read.cigarElements().stream().filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();
            baseTrimTotal += read.baseTrimCount();

            if(support.type().isSplitSupport())
            {
                softClipMatchTotal += support.junctionMatches();
                softClipMismatchTotal += support.junctionMismatches();
                refBaseMismatchTotal += support.referenceMismatches();

                if(support.junctionMatches() > maxExtBaseMatchCount)
                {
                    softClipSecondMaxLength = maxExtBaseMatchCount; // promote the second highest
                    maxExtBaseMatchCount = support.junctionMatches();
                }
                else if(support.junctionMatches() > softClipSecondMaxLength)
                {
                    softClipSecondMaxLength = support.junctionMatches();
                }
            }
        }

        sj.add(String.valueOf(softClipMatchTotal));
        sj.add(String.valueOf(softClipMismatchTotal));
        sj.add(String.valueOf(softClipSecondMaxLength));
        sj.add(String.valueOf(refBaseMismatchTotal));
        sj.add(String.valueOf(baseTrimTotal));

        int supportCount = assembly.supportCount();

        sj.add(statString(nmCountTotal, supportCount));
        sj.add(statString(indelLengthTotal, supportCount));
        sj.add(statString(baseQualTotal, supportCount));
        sj.add(statString(mapQualTotal, supportCount));
        sj.add(String.valueOf(assembly.mismatchReadCount()));
    }

    private static String statString(int count, double readCount)
    {
        double avgValue = count/readCount;
        return format("%d", round(avgValue));
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

        sj.add("SecondaryLinks");
    }

    protected static void addPhasingInfo(final JunctionAssembly assembly, final StringJoiner sj)
    {
        PhaseGroup phaseGroup = assembly.phaseGroup();

        if(phaseGroup != null)
        {
            sj.add(String.valueOf(phaseGroup.id()));
            sj.add(String.valueOf(phaseGroup.assemblyCount()));
        }
        else
        {
            sj.add("-1").add("0");
        }

        PhaseSet phaseSet = phaseGroup != null ? phaseGroup.findPhaseSet(assembly) : null;

        if(phaseSet != null)
        {
            sj.add(String.valueOf(phaseSet.id()));
            sj.add(String.valueOf(phaseSet.assemblies().size()));

            List<AssemblyLink> assemblyLinks = phaseSet.findAssemblyLinks(assembly);

            List<AssemblyLink> splitLinks = assemblyLinks.stream()
                    .filter(x -> x.type() == LinkType.SPLIT || x.type() == LinkType.INDEL).collect(Collectors.toList());

            List<AssemblyLink> facingLinks = assemblyLinks.stream().filter(x -> x.type() == LinkType.FACING).collect(Collectors.toList());
            sj.add(assemblyLinksStr(assembly, splitLinks));
            sj.add(assemblyLinksStr(assembly, facingLinks));

            if(!splitLinks.isEmpty())
            {
                AssemblyLink svLink = splitLinks.get(0);
                sj.add(svLink.svType().toString());
                sj.add(String.valueOf(svLink.length()));
                sj.add(svLink.insertedBases());
                sj.add(svLink.overlapBases());
            }
            else
            {
                sj.add("").add("0").add("").add("");
            }
        }
        else
        {
            sj.add("-1").add("0").add("").add("").add("").add("0").add("").add("");
        }

        List<AssemblyLink> secondarySplitLinks = phaseGroup != null ? phaseGroup.findSecondarySplitLinks(assembly) : Collections.emptyList();
        sj.add(assemblyLinksStr(assembly, secondarySplitLinks));
    }

    protected static String assemblyLinksStr(final JunctionAssembly assembly, final List<AssemblyLink> assemblyLinks)
    {
        if(assemblyLinks.isEmpty())
            return "";

        StringJoiner sj = new StringJoiner(ITEM_DELIM);

        for(AssemblyLink link : assemblyLinks)
        {
            JunctionAssembly otherAssembly = link.otherAssembly(assembly);
            sj.add(format("%s=%s", link.type(), otherAssembly.junction().coords()));
        }
        return sj.toString();
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
            remoteJunctMate += region.readTypeCounts()[REMOTE_READ_TYPE_JUNCTION_MATE];
            remoteJunctSupp += region.readTypeCounts()[REMOTE_READ_TYPE_JUNCTION_SUPP];
            remoteDiscordant += region.readTypeCounts()[REMOTE_READ_TYPE_DISCORDANT_READ];
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
