package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;
import static java.lang.String.valueOf;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.INDEL;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadUtils;

import org.jetbrains.annotations.Nullable;

public class AssemblyStats
{
    // read/fragment-type counts
    public int JuncSupps;
    public int JuncMateConcordant;
    public int JuncMateDiscordantRemote;
    public int JuncMateDiscordantRefSide;
    public int JuncMateRefSide;
    public int JuncMateUnmappedRemote;
    public int JuncMateUnmappedRefSide;
    public int IndelReads;

    public double ProximateJuncReadRatio;

    // read qualities
    public int IndelLengthTotal;
    public int BaseQualTotal;
    public int MapQualTotal;
    public int BaseTrimTotal;
    public int SoftClipMatchTotal;
    public int SoftClipMismatchTotal;
    public int RefBaseMismatchTotal;
    public int SoftClipSecondMaxLength;
    public int MaxExtBaseMatchCount;
    public int NonExtensionSupportReads;

    public int CandidateSupportCount;
    public int UnmappedReadCount;

    public int ReadCount;

    public final Map<String,Integer> RemoteSuppRegionFrequency;

    public AssemblyStats()
    {
        ReadCount = 0;
        JuncSupps = 0;
        JuncMateConcordant = 0;
        JuncMateDiscordantRemote = 0;
        JuncMateDiscordantRefSide = 0;
        JuncMateRefSide = 0;
        JuncMateUnmappedRemote = 0;
        JuncMateUnmappedRefSide = 0;
        IndelReads = 0;
        ProximateJuncReadRatio = 0;

        IndelLengthTotal = 0;
        BaseQualTotal = 0;
        MapQualTotal = 0;
        BaseTrimTotal = 0;
        SoftClipMatchTotal = 0;
        SoftClipMismatchTotal = 0;
        RefBaseMismatchTotal = 0;
        SoftClipSecondMaxLength = 0;
        MaxExtBaseMatchCount = 0;
        NonExtensionSupportReads = 0;
        CandidateSupportCount = 0;
        UnmappedReadCount = 0;
        RemoteSuppRegionFrequency = Maps.newHashMap();
    }

    public void addRead(final SupportRead supportRead, final Junction junction, @Nullable final Read read)
    {
        ++ReadCount;

        if(supportRead.type().isSplitSupport())
        {
            if(supportRead.type() == INDEL)
                ++IndelReads;
            else if(junction.indelBased() && supportRead.indelCoords() != null && supportRead.indelCoords().Length >= MIN_VARIANT_LENGTH)
                ++IndelReads;

            if(supportRead.isSupplementary())
            {
                ++JuncSupps;
            }
            else if(supportRead.isPairedRead())
            {
                // don't attempt to set mate info for supplementaries - these counts are for the primary reads (ie the fragment itself)
                boolean matePastJunction = (supportRead.orientation().isForward()) == junction.isForward();

                if(supportRead.isMateUnmapped())
                {
                    if(matePastJunction)
                        ++JuncMateUnmappedRemote;
                    else
                        ++JuncMateUnmappedRefSide;
                }
                else
                {
                    if(supportRead.isDiscordant())
                    {
                        if(matePastJunction)
                            ++JuncMateDiscordantRemote;
                        else
                            ++JuncMateDiscordantRefSide;
                    }
                    else
                    {
                        if((junction.isForward() && supportRead.mateAlignmentStart() > junction.Position)
                        || (!junction.isForward() && supportRead.mateAlignmentEnd() < junction.Position))
                        {
                            ++JuncMateConcordant;
                        }
                        else
                        {
                            ++JuncMateRefSide;
                        }
                    }
                }
            }

            SoftClipMatchTotal += supportRead.extensionBaseMatches();
            SoftClipMismatchTotal += supportRead.extensionBaseMismatches();
            RefBaseMismatchTotal += max(supportRead.referenceMismatches(), 0);

            if(supportRead.extensionBaseMatches() > MaxExtBaseMatchCount)
            {
                SoftClipSecondMaxLength = MaxExtBaseMatchCount; // promote the second highest
                MaxExtBaseMatchCount = supportRead.extensionBaseMatches();
            }
            else if(supportRead.extensionBaseMatches() > SoftClipSecondMaxLength)
            {
                SoftClipSecondMaxLength = supportRead.extensionBaseMatches();
            }

            // track frequency of remote supplementary locations
            if(junction.softClipBased())
                trackSuppRemoteRegions(supportRead, junction);
        }

        MapQualTotal += supportRead.mapQual();
        BaseTrimTotal += supportRead.trimCount();

        if(read != null)
        {
            BaseQualTotal += ReadUtils.avgBaseQuality(read);
            IndelLengthTotal += read.cigarElements().stream().filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();
        }
    }

    private void trackSuppRemoteRegions(final SupportRead read, final Junction junction)
    {
        if(read.supplementaryData() == null || read.supplementaryData().MapQuality < 30)
            return;

        if(read.junctionExtensionLength(junction.Orient) < ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
            return;

        // must soft-clip at the junction
        if(junction.isForward() && read.alignmentEnd() != junction.Position)
            return;
        else if(junction.isReverse() && read.alignmentStart() != junction.Position)
            return;

        String remoteLocationStr = format("%s_%d",
                read.supplementaryData().Chromosome, round(read.supplementaryData().Position / 1000.0));

        RemoteSuppRegionFrequency.merge(remoteLocationStr, 1, Integer::sum);
    }

    public double suppRemoteRegionRatio()
    {
        if(RemoteSuppRegionFrequency.isEmpty())
            return 0;

        if(RemoteSuppRegionFrequency.size() == 1)
            return 1.0;

        int total = 0;
        int maxRegion = 0;

        for(Integer count : RemoteSuppRegionFrequency.values())
        {
            total += count;
            maxRegion = max(maxRegion, count);
        }

        return total > 0 ? maxRegion / (double)total : 0;
    }

    public static void addReadTypeHeader(final StringJoiner sj)
    {
        sj.add("JuncSupps");
        sj.add("JuncIndels");
        sj.add("JuncMateConcord");
        sj.add("JuncMateRefSide");
        sj.add("JuncMateDiscordRemote");
        sj.add("JuncMateDiscordRefSide");
        sj.add("JuncMateUnmappedRemote");
        sj.add("JuncMateUnmappedRefSide");
    }

    public void addReadTypeCounts(final StringJoiner sj)
    {
        sj.add(valueOf(JuncSupps));
        sj.add(valueOf(IndelReads));

        sj.add(valueOf(JuncMateConcordant));
        sj.add(valueOf(JuncMateRefSide));
        sj.add(valueOf(JuncMateDiscordantRemote));
        sj.add(valueOf(JuncMateDiscordantRefSide));
        sj.add(valueOf(JuncMateUnmappedRemote));
        sj.add(valueOf(JuncMateUnmappedRefSide));
    }

    public static void addReadStatsHeader(final StringJoiner sj)
    {
        sj.add("SoftClipMatches");
        sj.add("SoftClipMismatches");
        sj.add("SoftClip2ndMaxLength");
        sj.add("RefBaseMismatches");
        sj.add("BaseTrimCount");
        sj.add("NonExtSupportReads");
        sj.add("ProxJuncReadRatio");
        sj.add("SuppRemoteLocationRatio");

        sj.add("AvgIndelLength");
        sj.add("AvgBaseQual");
        sj.add("AvgMapQual");
    }

    public void addReadStats(final StringJoiner sj)
    {
        sj.add(valueOf(SoftClipMatchTotal));
        sj.add(valueOf(SoftClipMismatchTotal));
        sj.add(valueOf(SoftClipSecondMaxLength));
        sj.add(valueOf(RefBaseMismatchTotal));
        sj.add(valueOf(BaseTrimTotal));
        sj.add(valueOf(NonExtensionSupportReads));
        sj.add(format("%.2f", ProximateJuncReadRatio));
        sj.add(format("%.2f", suppRemoteRegionRatio()));

        sj.add(statString(IndelLengthTotal, ReadCount));
        sj.add(statString(BaseQualTotal, ReadCount));
        sj.add(statString(MapQualTotal, ReadCount));
    }

    private static String statString(int count, double readCount)
    {
        double avgValue = count/readCount;
        return format("%d", round(avgValue));
    }
}
