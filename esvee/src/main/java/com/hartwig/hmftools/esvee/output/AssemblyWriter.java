package com.hartwig.hmftools.esvee.output;

import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_DISCORDANT_READ;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.common.RemoteRegion.REMOTE_READ_TYPE_JUNCTION_SUPP;
import static com.hartwig.hmftools.esvee.common.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.common.SupportType.INDEL;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION_MATE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.AssemblyLink;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.BaseMismatches;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.LinkType;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.PhaseSet;
import com.hartwig.hmftools.esvee.common.RefSideSoftClip;
import com.hartwig.hmftools.esvee.common.RemoteRegion;
import com.hartwig.hmftools.esvee.common.RepeatInfo;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadUtils;
import com.hartwig.hmftools.esvee.utils.TruthsetAnnotation;

public class AssemblyWriter
{
    private final SvConfig mConfig;

    private final BufferedWriter mWriter;
    private final TruthsetAnnotation mTruthsetAnnotation;

    // write info about assemblies
    public AssemblyWriter(final SvConfig config)
    {
        mConfig = config;
        mTruthsetAnnotation = new TruthsetAnnotation(mConfig.TruthsetFile);

        mWriter = initialiseWriter();
    }

    public void close() { closeBufferedWriter(mWriter);}

    private BufferedWriter initialiseWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.ASSEMBLIES))
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.outputFilename(WriteType.ASSEMBLIES));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("Id");
            sj.add("Chromosome").add("JunctionPosition").add("JunctionOrientation");

            sj.add("ExtBaseLength").add("RefBasePosition").add("RefBaseLength");

            sj.add("SplitFrags").add("RefSplitFrags");
            sj.add("DiscFrags").add("RefDiscFrags");

            sj.add("IndelReads");
            sj.add("JuncMates");
            sj.add("JuncMateUnmapped");

            sj.add("SoftClipMismatches");
            sj.add("RefBaseMismatches");

            sj.add("AvgNmCount");
            sj.add("AvgIndelLength");
            sj.add("AvgBaseQual");
            sj.add("AvgMapQual");
            sj.add("AvgTrimCount");
            sj.add("InitialReadId");

            sj.add("RepeatInfo");

            sj.add("RefSideSoftClips");

            sj.add("RemoteRegionCount");
            sj.add("RemoteRegionInfo");
            sj.add("RemoteRegionJuncMate");
            sj.add("RemoteRegionJuncSupps");
            sj.add("RemoteRegionDiscordant");

            sj.add("JuncUnlinkedMates");
            sj.add("JuncUnlinkedSupps");
            sj.add("DiscUnlinkedMates");
            sj.add("InitRefBaseCandidates");

            sj.add("Outcome");

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

            sj.add("JunctionSequence");
            sj.add("RefBaseSequence");

            sj.add("MergedAssemblies");
            sj.add("BranchedAssemblyIds");

            if(mTruthsetAnnotation.enabled())
                sj.add(TruthsetAnnotation.tsvHeader());

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise assembly writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeAssembly(final JunctionAssembly assembly)
    {
        if(mWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(String.valueOf(assembly.id()));
            sj.add(assembly.junction().Chromosome);
            sj.add(String.valueOf(assembly.junction().Position));
            sj.add(String.valueOf(assembly.junction().Orientation));

            sj.add(String.valueOf(assembly.extensionLength()));
            sj.add(String.valueOf(assembly.refBasePosition()));
            sj.add(String.valueOf(assembly.refBaseLength()));

            int juncFrags = 0;
            int refSampleJuncFrags = 0;
            int discFrags = 0;
            int refSampleDiscFrags = 0;
            int indelReads = 0;
            int juncMateUnmapped = 0;
            int juncMateCount = 0;
            int refBaseMismatches = 0;
            int softClipBaseMismatches = 0;

            int juncUnlinkedMates = 0;
            int juncUnlinkedSupps = 0;
            int discUnlinkedMates = 0;

            Set<String> uniqueReadIds = Sets.newHashSet();

            for(AssemblySupport support : assembly.support())
            {
                boolean isReference = support.read().isReference();
                Read read = support.read();

                if(support.type() == JUNCTION || support.type() == INDEL)
                {
                    if(uniqueReadIds.contains(read.getName()))
                        continue;

                    uniqueReadIds.add(read.getName()); // since fragments are being counted

                    ++juncFrags;

                    if(support.type() == INDEL)
                        ++indelReads;

                    if(isReference)
                        ++refSampleJuncFrags;

                    softClipBaseMismatches += support.junctionMismatches();
                    refBaseMismatches += support.referenceMismatches();

                    if(read.isMateUnmapped())
                        ++juncMateUnmapped;

                    if(read.isMateMapped() && !read.isSupplementary() && !read.hasMateSet())
                        ++juncUnlinkedMates;

                    if(read.hasSupplementary() && read.supplementaryRead() == null)
                        ++juncUnlinkedSupps;
                }
                else
                {
                    refBaseMismatches += support.mismatchCount();

                    if(support.type() == DISCORDANT)
                    {
                        ++discFrags;

                        if(isReference)
                            ++refSampleDiscFrags;

                        if(!read.isSupplementary() && !read.hasMateSet())
                            ++discUnlinkedMates;
                    }
                    else if(support.type() == JUNCTION_MATE)
                    {
                        ++juncMateCount;
                    }
                }
            }

            sj.add(String.valueOf(juncFrags));
            sj.add(String.valueOf(refSampleJuncFrags));
            sj.add(String.valueOf(discFrags));
            sj.add(String.valueOf(refSampleDiscFrags));
            sj.add(String.valueOf(indelReads));
            sj.add(String.valueOf(juncMateCount));
            sj.add(String.valueOf(juncMateUnmapped));

            sj.add(String.valueOf(softClipBaseMismatches));
            sj.add(String.valueOf(refBaseMismatches));

            // ref sequence stats purely for analysis
            ReadStats readStats = buildReadStats(assembly.support());
            sj.add(statString(readStats.NmCountTotal, assembly.supportCount()));
            sj.add(statString(readStats.IndelLengthTotal, assembly.supportCount()));
            sj.add(statString(readStats.BaseQualTotal, assembly.supportCount()));
            sj.add(statString(readStats.MapQualTotal, assembly.supportCount()));
            sj.add(statString(readStats.BaseTrimCountTotal, assembly.supportCount()));

            sj.add(assembly.initialRead() != null ? assembly.initialRead().getName() : "NONE"); // shouldn't occur

            sj.add(repeatsInfoStr(assembly.repeatInfo()));

            sj.add(refSideSoftClipsStr(assembly.refSideSoftClips()));

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

            sj.add(String.valueOf(juncUnlinkedMates));
            sj.add(String.valueOf(juncUnlinkedSupps));
            sj.add(String.valueOf(discUnlinkedMates));
            sj.add(String.valueOf(assembly.candidateSupport().size()));

            sj.add(String.valueOf(assembly.outcome()));

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

            if(assembly.hasUnsetBases())
            {
                sj.add("UNSET_BASES");
                sj.add("UNSET_BASES");
            }
            else
            {
                sj.add(assembly.formJunctionSequence());

                int refBaseLength = mConfig.AssemblyRefBaseWriteMax == 0 ? assembly.refBaseLength() : mConfig.AssemblyRefBaseWriteMax;
                sj.add(assembly.formRefBaseSequence(refBaseLength)); // long enough to show most short TIs
            }

            sj.add(String.valueOf(assembly.mergedAssemblyCount()));

            String branchedAssemblyIds = assembly.branchedAssemblies().stream()
                    .map(x -> String.valueOf(x.id())).collect(Collectors.joining(";"));
            sj.add(branchedAssemblyIds);

            if(mTruthsetAnnotation.enabled())
                sj.add(mTruthsetAnnotation.findTruthsetAnnotation(assembly));

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write assembly: {}", e.toString());
        }
    }

    private static String statString(int count, double readCount)
    {
        double avgValue = count/readCount;
        return format("%d", round(avgValue));
    }

    private class ReadStats
    {
        public int NmCountTotal;
        public int IndelLengthTotal;
        public int BaseQualTotal;
        public int MapQualTotal;
        public int BaseTrimCountTotal;

        public ReadStats()
        {
            NmCountTotal = 0;
            IndelLengthTotal = 0;
            BaseQualTotal = 0;
            MapQualTotal = 0;
            BaseTrimCountTotal = 0;
        }
    }

    private ReadStats buildReadStats(final List<AssemblySupport> supportReads)
    {
        ReadStats readStats = new ReadStats();

        for(AssemblySupport support : supportReads)
        {
            Read read = support.read();
            readStats.NmCountTotal += read.numberOfEvents();
            readStats.MapQualTotal += read.mappingQuality();
            readStats.BaseQualTotal += ReadUtils.avgBaseQuality(read);
            readStats.IndelLengthTotal += read.cigarElements().stream().filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();
            readStats.BaseTrimCountTotal += read.baseTrimCount();
        }

        return readStats;
    }

    private static String refSideSoftClipsStr(final List<RefSideSoftClip> refSideSoftClips)
    {
        if(refSideSoftClips.isEmpty())
            return "";

        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        refSideSoftClips.forEach(x -> sj.add(format("%d:%d=%d", x.Position, x.maxLength(), x.readCount())));
        return sj.toString();
    }

    private static String assemblyLinksStr(final JunctionAssembly assembly, final List<AssemblyLink> assemblyLinks)
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

    private static String repeatsInfoStr(final List<RepeatInfo> repeats)
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

}
