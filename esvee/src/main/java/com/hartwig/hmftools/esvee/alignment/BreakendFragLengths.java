package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

public class BreakendFragLengths
{
    public BreakendFragLengths() { }

    public void calcAssemblyFragmentLengths(final List<List<AssemblyAlignment>> assemblyAlignmentGroups)
    {
        assemblyAlignmentGroups.forEach(x -> calcInferredFragmentLengths(x));
    }

    private void calcInferredFragmentLengths(final List<AssemblyAlignment> assemblyAlignments)
    {
        Map<String, BreakendFragmentData> fragmentSupportMap = Maps.newHashMap();

        for(AssemblyAlignment assemblyAlignment : assemblyAlignments)
        {
            for(Breakend breakend : assemblyAlignment.breakends())
            {
                for(JunctionAssembly assembly : assemblyAlignment.assemblies())
                {
                    for(SupportRead read : assembly.support())
                    {
                        if(read.isSupplementary())
                            continue;

                        if(!breakend.Chromosome.equals(read.chromosome()))
                            continue;

                        boolean isSplitFragment = false;
                        boolean isDiscFragment = false;

                        if(breakend.readSpansJunction(read, false))
                        {
                            isSplitFragment = true;
                        }
                        else if(breakend.isRelatedDiscordantRead(read.alignmentStart(), read.alignmentEnd(), read.orientation()))
                        {
                            isDiscFragment = true;
                        }

                        if(!isSplitFragment && !isDiscFragment)
                            continue;

                        BreakendFragmentData fragmentSupport = fragmentSupportMap.get(read.id());

                        if(fragmentSupport == null)
                        {
                            fragmentSupport = new BreakendFragmentData(read, breakend);
                            fragmentSupportMap.put(read.id(), fragmentSupport);
                        }
                        else
                        {
                            fragmentSupport.add(read, breakend);
                        }
                    }
                }
            }
        }

        Map<Breakend,LengthData> breakendFragmentLengths = Maps.newHashMap();

        for(Map.Entry<String,BreakendFragmentData> entry : fragmentSupportMap.entrySet())
        {
            BreakendFragmentData breakendFragmentData = entry.getValue();

            int inferredFragmentLength = calcInferredFragmentLength(breakendFragmentData);
            if(inferredFragmentLength != INVALID_FRAGMENT_LENGTH)
            {
                breakendFragmentData.Reads.forEach(x -> x.setInferredFragmentLength(inferredFragmentLength));

                for(Breakend breakend : breakendFragmentData.Breakends)
                {
                    LengthData lengthData = breakendFragmentLengths.get(breakend);

                    if(lengthData == null)
                    {
                        lengthData = new LengthData();
                        breakendFragmentLengths.put(breakend, lengthData);
                    }

                    ++lengthData.Count;
                    lengthData.TotalLength += inferredFragmentLength;
                }
            }
        }

        for(Map.Entry<Breakend,LengthData> entry : breakendFragmentLengths.entrySet())
        {
            Breakend breakend = entry.getKey();
            LengthData lengthData = entry.getValue();
            breakend.setAverageInferredFragmentLength(lengthData.averageLength());
        }
    }

    private class LengthData
    {
        public int Count;
        public int TotalLength;

        public LengthData()
        {
            Count = 0;
            TotalLength = 0;
        }

        public int averageLength() { return (int)round(TotalLength / (double)Count); }
    }

    public static final int INVALID_FRAGMENT_LENGTH = -1;

    private static Breakend findRelatedBreakend(final SupportRead read, final List<Breakend> breakends)
    {
        for(Breakend breakend : breakends)
        {
            if(breakend.readSpansJunction(read, false))
                return breakend;

            if(breakend.isRelatedDiscordantRead(read.alignmentStart(), read.alignmentEnd(), read.orientation()))
                return breakend;
        }

        return null;
    }

    public static int calcInferredFragmentLength(final BreakendFragmentData breakendFragmentData)
    {
        if(breakendFragmentData.Reads.size() != 2 || breakendFragmentData.Breakends.isEmpty())
            return INVALID_FRAGMENT_LENGTH;

        SupportRead read1 = breakendFragmentData.Reads.get(0);
        SupportRead read2 = breakendFragmentData.Reads.get(1);

        Breakend breakend1 = findRelatedBreakend(read1, breakendFragmentData.Breakends);
        Breakend breakend2 = findRelatedBreakend(read2, breakendFragmentData.Breakends);

        if(breakend1 == null || breakend2 == null)
            return INVALID_FRAGMENT_LENGTH;

        if(breakend1 == breakend2)
        {
            return max(read1.unclippedEnd(), read2.unclippedEnd()) - min(read1.unclippedStart(), read2.unclippedStart());
        }

        int junctionDistance = breakend1.Orient.isForward() ?
                breakend1.Position - read1.unclippedStart() : read1.unclippedEnd() - breakend1.Position;

        if(breakend1.otherBreakend() == breakend2)
        {
            junctionDistance += breakend2.Orient.isForward() ?
                    breakend2.Position - read2.unclippedStart() : read2.unclippedEnd() - breakend2.Position;

            // immediately across a junction, so take the distance from each
            return junctionDistance;
        }

        // factor in chained breakend links
        Breakend prevBreakend = breakend1;
        Breakend nextBreakend;
        boolean nextIsFacing = false;

        if(read1.orientation() == breakend1.Orient)
        {
            nextBreakend = breakend1.otherBreakend();
            nextIsFacing = false;
        }
        else
        {
            nextBreakend = !breakend1.facingBreakends().isEmpty() ? breakend1.facingBreakends().get(0) : null;
            nextIsFacing = true;
        }

        Set<Breakend> processedBreakends = Sets.newHashSet(breakend1);

        while(nextBreakend != null)
        {
            if(processedBreakends.contains(nextBreakend))
            {
                SV_LOGGER.debug("inferred frag length breakend({}) re-encountered for read({})", nextBreakend, read1);
                break;
            }

            processedBreakends.add(nextBreakend);

            if(nextBreakend == breakend2)
            {
                if(breakend2.Orient.isForward())
                {
                    if(nextIsFacing)
                        junctionDistance += max(breakend2.Position, read2.unclippedEnd()) - prevBreakend.Position;
                    else
                        junctionDistance += breakend2.Position - read2.unclippedStart();
                }
                else
                {
                    if(nextIsFacing)
                        junctionDistance += prevBreakend.Position - min(breakend2.Position, read2.unclippedStart());
                    else
                        junctionDistance += read2.unclippedEnd() - breakend2.Position;
                }

                break;
            }

            if(nextIsFacing)
            {
                junctionDistance += abs(prevBreakend.Position - nextBreakend.Position);
                prevBreakend = nextBreakend;
                nextBreakend = nextBreakend.otherBreakend();
                nextIsFacing = false;
            }
            else
            {
                prevBreakend = nextBreakend;
                nextBreakend = nextBreakend.facingBreakends().stream().filter(x -> !processedBreakends.contains(x)).findFirst().orElse(null);
                nextIsFacing = true;
            }
        }

        return junctionDistance;
    }

    private class BreakendFragmentData
    {
        public final List<Breakend> Breakends;
        public final List<SupportRead> Reads;

        public BreakendFragmentData(final SupportRead read, final Breakend breakend)
        {
            Breakends = Lists.newArrayList(breakend);
            Reads = Lists.newArrayList(read);
        }

        public void add(final SupportRead read, final Breakend breakend)
        {
            if(Reads.stream().noneMatch(x -> x.flags() == read.flags()))
                Reads.add(read);

            if(!Breakends.contains(breakend))
                Breakends.add(breakend);
        }

        public String toString() { return format("breakends(%d) reads(%d)", Breakends.size(), Reads.size()); }
    }
}
