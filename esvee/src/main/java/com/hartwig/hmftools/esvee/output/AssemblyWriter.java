package com.hartwig.hmftools.esvee.output;

import static java.lang.Math.max;
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
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION_MATE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.BaseMismatches;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.RefBaseAssembly;
import com.hartwig.hmftools.esvee.common.RefSideSoftClip;
import com.hartwig.hmftools.esvee.common.RemoteRegion;
import com.hartwig.hmftools.esvee.common.RepeatInfo;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadUtils;

public class AssemblyWriter
{
    private final SvConfig mConfig;

    private final BufferedWriter mWriter;

    // write info about assemblies
    public AssemblyWriter(final SvConfig config)
    {
        mConfig = config;

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
            sj.add("Chromosome");
            sj.add("JunctionPosition");
            sj.add("JunctionOrientation");
            sj.add("SoftClipLength");
            sj.add("RefBasePosition");
            sj.add("RefBaseLength");
            sj.add("SupportCount");
            sj.add("RefSupportCount");
            sj.add("SoftClipMismatches");
            sj.add("RefBaseMismatches");
            sj.add("RefBaseDominantMismatches");

            sj.add("AvgNmCount");
            sj.add("AvgIndelLength");
            sj.add("AvgBaseQual");
            sj.add("AvgMapQual");
            sj.add("InitialReadId");

            sj.add("RepeatInfo");

            sj.add("RefSideSoftClips");

            sj.add("RefExtDistance");
            sj.add("RefExtJuncMates");
            sj.add("RefExtDiscReads");
            sj.add("RefExtMismatches");

            sj.add("UnmappedJuncReads");

            sj.add("RemoteRegionCount");
            sj.add("RemoteRegionJuncMate");
            sj.add("RemoteRegionJuncSupps");
            sj.add("RemoteRegionDiscordant");
            sj.add("RemoteRegionInfo");

            sj.add("PhaseGroupId");
            sj.add("PhaseGroupCount");
            sj.add("MergedAssemblies");
            sj.add("BranchedAssemblyIds");

            sj.add("JunctionSequence");
            sj.add("RefBaseSequence");

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

            if(assembly.refBaseAssembly() != null)
            {
                sj.add(String.valueOf(assembly.refBaseAssembly().extensionRefPosition()));
                sj.add(String.valueOf(assembly.refBaseAssembly().baseLength())); // since includes the junction base itself

            }
            else
            {
                sj.add(String.valueOf(assembly.isForwardJunction() ? assembly.minAlignedPosition() : assembly.maxAlignedPosition()));
                sj.add(String.valueOf(assembly.refBaseLength()));
            }

            int refBaseMismatches = 0;
            int softClipBaseMismatches = 0;

            for(AssemblySupport support : assembly.support())
            {
                refBaseMismatches += support.referenceMismatches();
                softClipBaseMismatches += support.junctionMismatches();
            }

            // where the mismatches on a ref base exceeds 50% of the junction read count, suggesting the wrong base was used or there are
            // valid alternatives
            int refBaseDominantMismatches = 0;

            if(assembly.mismatches().hasMismatches())
            {
                for(BaseMismatches baseMismatch : assembly.mismatches().indexedBaseMismatches().values())
                {
                    if(baseMismatch.mismatchReadTotal() >= assembly.supportCount() * 0.4)
                        ++refBaseDominantMismatches;
                }
            }

            sj.add(String.valueOf(assembly.supportCount()));
            long referenceSupportCount = assembly.support().stream().filter(x -> x.read().isReference()).count();
            sj.add(String.valueOf(referenceSupportCount));
            sj.add(String.valueOf(softClipBaseMismatches));
            sj.add(String.valueOf(refBaseMismatches));
            sj.add(String.valueOf(refBaseDominantMismatches));

            // ref sequence stats purely for analysis
            ReadStats readStats = buildReadStats(assembly.support(), assembly.junction().isForward());
            sj.add(statString(readStats.NmCountTotal, assembly.supportCount()));
            sj.add(statString(readStats.IndelLengthTotal, assembly.supportCount()));
            sj.add(statString(readStats.BaseQualTotal, assembly.supportCount()));
            sj.add(statString(readStats.MapQualTotal, assembly.supportCount()));

            sj.add(assembly.initialRead().getName());

            sj.add(repeatsInfoStr(assembly.repeatInfo()));

            sj.add(refSideSoftClipsStr(assembly.refSideSoftClips()));

            RefBaseAssembly refBaseAssembly = assembly.refBaseAssembly();

            if(refBaseAssembly != null)
            {
                sj.add(String.valueOf(refBaseAssembly.nonJunctionReadExtension()));

                int juncMateCount = (int)refBaseAssembly.support().stream().filter(x -> x.type() == JUNCTION_MATE).count();
                int discordantReadCount = (int)refBaseAssembly.support().stream().filter(x -> x.type() == DISCORDANT).count();
                sj.add(String.valueOf(juncMateCount));
                sj.add(String.valueOf(discordantReadCount));
                sj.add(String.valueOf(refBaseAssembly.mismatches().readBaseCount()));
            }
            else
            {
                sj.add("0").add("0").add("0").add("0");
            }

            int unmappedJuncReads = (int)assembly.support().stream().filter(x -> x.type() == JUNCTION && x.read().isMateUnmapped()).count();
            sj.add(String.valueOf(unmappedJuncReads));

            sj.add(String.valueOf(assembly.remoteRegions().size()));

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

            sj.add(remoteRegionInfoStr(assembly.remoteRegions()));

            if(assembly.primaryPhaseGroup() != null)
            {
                sj.add(String.valueOf(assembly.primaryPhaseGroup().id()));
                sj.add(String.valueOf(assembly.primaryPhaseGroup().assemblyCount()));
            }
            else
            {
                sj.add("-1").add("0");
            }

            sj.add(String.valueOf(assembly.mergedAssemblyCount()));

            String branchedAssemblyIds = assembly.branchedAssemblies().stream()
                    .map(x -> String.valueOf(x.id())).collect(Collectors.joining(";"));
            sj.add(branchedAssemblyIds);

            sj.add(assembly.formJunctionSequence());
            sj.add(assembly.formRefBaseSequence());

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

        public ReadStats()
        {
            NmCountTotal = 0;
            IndelLengthTotal = 0;
            BaseQualTotal = 0;
            MapQualTotal = 0;
        }
    }

    private ReadStats buildReadStats(final List<AssemblySupport> supportReads, boolean isForwardJunction)
    {
        ReadStats readStats = new ReadStats();

        for(AssemblySupport support : supportReads)
        {
            Read read = support.read();
            readStats.NmCountTotal += read.numberOfEvents();
            readStats.MapQualTotal += read.mappingQuality();
            readStats.BaseQualTotal += ReadUtils.avgBaseQuality(read);
            readStats.IndelLengthTotal += read.cigarElements().stream().filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();
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
