package com.hartwig.hmftools.esvee.output;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

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

import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.WriteType;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.BaseMismatches;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.RefBaseAssembly;
import com.hartwig.hmftools.esvee.common.RemoteRegion;
import com.hartwig.hmftools.esvee.common.RepeatInfo;
import com.hartwig.hmftools.esvee.common.SupportType;
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

            sj.add("Chromosome");
            sj.add("JunctionPosition");
            sj.add("JunctionOrientation");
            sj.add("RangeStart");
            sj.add("RangeEnd");
            sj.add("Length");
            sj.add("SupportCount");
            sj.add("RefSupportCount");
            sj.add("SoftClipMismatches");
            sj.add("RefBaseMismatches");
            sj.add("RefBaseDominantMismatches");
            sj.add("JunctionSequence");

            sj.add("AvgNmCount");
            sj.add("AvgIndelLength");
            sj.add("AvgRefSideSoftClip");
            sj.add("AvgBaseQual");
            sj.add("AvgMapQual");
            sj.add("InitialReadId");

            sj.add("RepeatInfo");
            sj.add("MergedAssemblies");

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

            sj.add(assembly.junction().Chromosome);
            sj.add(String.valueOf(assembly.junction().Position));
            sj.add(String.valueOf(assembly.junction().Orientation));
            sj.add(String.valueOf(assembly.minAlignedPosition()));
            sj.add(String.valueOf(assembly.maxAlignedPosition()));
            sj.add(String.valueOf(assembly.baseLength()));

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

            sj.add(assembly.formSequence(5));

            // ref sequence stats purely for analysis
            ReadStats readStats = buildReadStats(assembly.support(), assembly.junction().isForward());
            sj.add(statString(readStats.NmCountTotal, assembly.supportCount()));
            sj.add(statString(readStats.IndelLengthTotal, assembly.supportCount()));
            sj.add(statString(readStats.RefSideSoftClipLengthTotal, assembly.supportCount()));
            sj.add(statString(readStats.BaseQualTotal, assembly.supportCount()));
            sj.add(statString(readStats.MapQualTotal, assembly.supportCount()));

            sj.add(assembly.initialRead().getName());

            sj.add(repeatsInfoStr(assembly.repeatInfo()));

            sj.add(String.valueOf(assembly.mergedAssemblyCount()));

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

            int unmappedJuncReads = (int)assembly.support().stream().filter(x -> x.type() == JUNCTION && !x.read().isMateMapped()).count();
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
        public int RefSideSoftClipLengthTotal;
        public int BaseQualTotal;
        public int MapQualTotal;

        public ReadStats()
        {
            NmCountTotal = 0;
            IndelLengthTotal = 0;
            RefSideSoftClipLengthTotal = 0;
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

            if(isForwardJunction)
                readStats.RefSideSoftClipLengthTotal += read.alignmentStart() - read.unclippedStart();
            else
                readStats.RefSideSoftClipLengthTotal += read.unclippedEnd() - read.alignmentEnd();
        }

        return readStats;
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
            sj.add(format("%s=%d", regions.get(i), regions.get(i).readCount()));
        }

        return sj.toString();
    }

}
