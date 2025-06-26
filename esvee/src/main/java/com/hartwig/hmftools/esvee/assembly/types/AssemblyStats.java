package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.types.SupportType.INDEL;

import java.util.StringJoiner;

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

    public int CandidateSupportCount;
    public int UnmappedReadCount;

    public int ReadCount;

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

        IndelLengthTotal = 0;
        BaseQualTotal = 0;
        MapQualTotal = 0;
        BaseTrimTotal = 0;
        SoftClipMatchTotal = 0;
        SoftClipMismatchTotal = 0;
        RefBaseMismatchTotal = 0;
        SoftClipSecondMaxLength = 0;
        MaxExtBaseMatchCount = 0;
        CandidateSupportCount = 0;
        UnmappedReadCount = 0;
    }

    public void addRead(final SupportRead supportRead, final Junction junction, @Nullable final Read read)
    {
        ++ReadCount;

        if(supportRead.type().isSplitSupport())
        {
            if(supportRead.type() == INDEL || supportRead.hasIndel())
                ++IndelReads;

            if(supportRead.isSupplementary())
            {
                ++JuncSupps;
            }
            else
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
        }

        MapQualTotal += supportRead.mapQual();
        BaseTrimTotal += supportRead.trimCount();

        if(read != null)
        {
            BaseQualTotal += ReadUtils.avgBaseQuality(read);
            IndelLengthTotal += read.cigarElements().stream().filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();
        }

        if(supportRead.type().isSplitSupport())
        {
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
        }
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
        sj.add(String.valueOf(JuncSupps));
        sj.add(String.valueOf(IndelReads));

        sj.add(String.valueOf(JuncMateConcordant));
        sj.add(String.valueOf(JuncMateRefSide));
        sj.add(String.valueOf(JuncMateDiscordantRemote));
        sj.add(String.valueOf(JuncMateDiscordantRefSide));
        sj.add(String.valueOf(JuncMateUnmappedRemote));
        sj.add(String.valueOf(JuncMateUnmappedRefSide));
    }

    public static void addReadStatsHeader(final StringJoiner sj)
    {
        sj.add("SoftClipMatches");
        sj.add("SoftClipMismatches");
        sj.add("SoftClip2ndMaxLength");
        sj.add("RefBaseMismatches");
        sj.add("BaseTrimCount");

        sj.add("AvgIndelLength");
        sj.add("AvgBaseQual");
        sj.add("AvgMapQual");
    }

    public void addReadStats(final StringJoiner sj)
    {
        sj.add(String.valueOf(SoftClipMatchTotal));
        sj.add(String.valueOf(SoftClipMismatchTotal));
        sj.add(String.valueOf(SoftClipSecondMaxLength));
        sj.add(String.valueOf(RefBaseMismatchTotal));
        sj.add(String.valueOf(BaseTrimTotal));

        sj.add(statString(IndelLengthTotal, ReadCount));
        sj.add(statString(BaseQualTotal, ReadCount));
        sj.add(statString(MapQualTotal, ReadCount));
    }

    public double avgMapQual() { return ReadCount > 0 ? MapQualTotal / (double)ReadCount : 0; }

    private static String statString(int count, double readCount)
    {
        double avgValue = count/readCount;
        return format("%d", round(avgValue));
    }
}
