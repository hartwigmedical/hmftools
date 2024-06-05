package com.hartwig.hmftools.esvee.assembly.types;

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
    public int NmCountTotal;
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

        NmCountTotal = 0;
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
                        // check if the mate

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

        NmCountTotal += supportRead.numOfEvents();
        MapQualTotal += supportRead.mapQual();
        BaseTrimTotal += supportRead.trimCount();

        if(read != null)
        {
            BaseQualTotal += ReadUtils.avgBaseQuality(read);
            IndelLengthTotal += read.cigarElements().stream().filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();
        }

        if(supportRead.type().isSplitSupport())
        {
            SoftClipMatchTotal += supportRead.junctionMatches();
            SoftClipMismatchTotal += supportRead.junctionMismatches();
            RefBaseMismatchTotal += supportRead.referenceMismatches();

            if(supportRead.junctionMatches() > MaxExtBaseMatchCount)
            {
                SoftClipSecondMaxLength = MaxExtBaseMatchCount; // promote the second highest
                MaxExtBaseMatchCount = supportRead.junctionMatches();
            }
            else if(supportRead.junctionMatches() > SoftClipSecondMaxLength)
            {
                SoftClipSecondMaxLength = supportRead.junctionMatches();
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

        sj.add("AvgNmCount");
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

        sj.add(statString(NmCountTotal, ReadCount));
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
