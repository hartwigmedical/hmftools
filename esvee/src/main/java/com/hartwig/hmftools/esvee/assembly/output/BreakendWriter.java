package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENTATION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.CommonUtils.withinLineProximity;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment;
import com.hartwig.hmftools.esvee.assembly.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.alignment.Breakend;
import com.hartwig.hmftools.esvee.assembly.alignment.BreakendSegment;
import com.hartwig.hmftools.esvee.assembly.types.InsertionType;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.WriteType;

public class BreakendWriter
{
    private final AssemblyConfig mConfig;

    private final BufferedWriter mWriter;

    public BreakendWriter(final AssemblyConfig config)
    {
        mConfig = config;

        mWriter = initialiseWriter();
    }

    public void close() { closeBufferedWriter(mWriter);}

    public static String FLD_BREAKEND_INS_SEQ = "InsertedBases";
    public static String FLD_BREAKEND_MATE_CHR = "MateChr";
    public static String FLD_BREAKEND_MATE_POSITION = "MatePos";
    public static String FLD_BREAKEND_MATE_ORIENT = "MateOrient";
    public static String FLD_SV_TYPE = "Type";

    private BufferedWriter initialiseWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.BREAKEND))
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.outputFilename(WriteType.BREAKEND));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("Id");
            sj.add("PhaseGroupId");
            sj.add("PhaseSetId");
            sj.add("AssemblyId");
            sj.add("MateId");
            sj.add("AssemblyInfo");

            sj.add(FLD_SV_TYPE).add(FLD_CHROMOSOME).add(FLD_POSITION).add(FLD_ORIENTATION);

            sj.add(FLD_BREAKEND_MATE_CHR).add(FLD_BREAKEND_MATE_POSITION).add(FLD_BREAKEND_MATE_ORIENT).add("Length");
            sj.add(FLD_BREAKEND_INS_SEQ).add("Homology").add("ConfidenceInterval").add("InexactOffset");

            sj.add("Qual");
            sj.add("SplitFragments");
            sj.add("RefSplitFragments");
            sj.add("DiscFragments");
            sj.add("RefDiscFragments");
            sj.add("ForwardReads");
            sj.add("ReverseReads");

            sj.add("SequenceLength");
            sj.add("SegmentCount");
            sj.add("SegmentIndex");
            sj.add("SequenceIndex");
            sj.add("AlignedBases");
            sj.add("MapQual");
            sj.add("Score");
            sj.add("AdjAlignedBases");
            sj.add("AvgFragmentLength");
            sj.add("IncompleteFragments");
            sj.add("BreakendQual");

            sj.add("FacingBreakendIds");

            sj.add("AltAlignments");
            sj.add("InsertionType");
            sj.add("UniqueFragPos");
            sj.add("ClosestAssembly");
            sj.add("NonPrimaryFragments");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise breakend writer: {}", e.toString());
            return null;
        }
    }

    public void writeBreakends(final AssemblyAlignment assemblyAlignment)
    {
        if(mWriter == null)
            return;

        try
        {
            String assemblyInfo = assemblyAlignment.info();

            Map<Breakend,String> closestAssemblyMap = Maps.newHashMap();

            for(Breakend breakend : assemblyAlignment.breakends())
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);

                sj.add(String.valueOf(breakend.id()));
                sj.add(String.valueOf(assemblyAlignment.assemblies().get(0).phaseGroup().id()));
                sj.add(String.valueOf(assemblyAlignment.phaseSet() != null ? assemblyAlignment.phaseSet().id() : -1));
                sj.add(String.valueOf(assemblyAlignment.id()));
                sj.add(!breakend.isSingle() ? String.valueOf(breakend.otherBreakend().id()) : "");
                sj.add(assemblyInfo);
                sj.add(String.valueOf(breakend.svType()));
                sj.add(breakend.Chromosome);
                sj.add(String.valueOf(breakend.Position));
                sj.add(String.valueOf(breakend.Orient));

                if(breakend.otherBreakend() != null)
                {
                    sj.add(breakend.otherBreakend().Chromosome);
                    sj.add(String.valueOf(breakend.otherBreakend().Position));
                    sj.add(String.valueOf(breakend.otherBreakend().Orient));
                    sj.add(String.valueOf(breakend.svLength()));

                }
                else
                {
                    sj.add("").add("").add("").add("0");
                }

                sj.add(breakend.InsertedBases);

                if(breakend.Homology.exists())
                {
                    sj.add(breakend.Homology.Homology);
                    sj.add(format("%d,%d", breakend.Homology.ExactStart, breakend.Homology.ExactEnd));
                    sj.add(format("%d,%d", breakend.Homology.InexactStart, breakend.Homology.InexactEnd));
                }
                else
                {
                    sj.add("").add("0,0").add("0,0");
                }

                sj.add(String.valueOf(breakend.calcSvQual()));

                int tumorCount = mConfig.TumorIds.size();

                int splitFrags = 0;
                int refSplitFrags = 0;
                int discFrags = 0;
                int refDiscFrags = 0;
                int forwardReads = 0;
                int reverseReads = 0;

                for(int i = 0; i < breakend.sampleSupport().size(); ++i)
                {
                    splitFrags += breakend.sampleSupport().get(i).SplitFragments;
                    discFrags += breakend.sampleSupport().get(i).DiscordantFragments;

                    if(i >= tumorCount)
                    {
                        refSplitFrags += breakend.sampleSupport().get(i).SplitFragments;
                        refDiscFrags += breakend.sampleSupport().get(i).DiscordantFragments;
                    }

                    forwardReads += breakend.sampleSupport().get(i).ForwardReads;
                    reverseReads += breakend.sampleSupport().get(i).ReverseReads;
                }

                sj.add(String.valueOf(splitFrags));
                sj.add(String.valueOf(refSplitFrags));
                sj.add(String.valueOf(discFrags));
                sj.add(String.valueOf(refDiscFrags));
                sj.add(String.valueOf(forwardReads));
                sj.add(String.valueOf(reverseReads));

                sj.add(String.valueOf(assemblyAlignment.fullSequenceLength()));

                // for now just the first segment - no showing branching or duplicates
                BreakendSegment segment = breakend.segments().get(0);
                sj.add(String.valueOf(breakend.segments().size()));
                sj.add(String.valueOf(segment.Index));
                sj.add(String.valueOf(segment.SequenceIndex));
                sj.add(String.valueOf(segment.Alignment.alignedBases()));
                sj.add(String.valueOf(segment.Alignment.mapQual()));
                sj.add(String.valueOf(segment.Alignment.score()));
                sj.add(String.valueOf(segment.Alignment.adjustedAlignment()));
                sj.add(String.valueOf(breakend.averageFragmentLength()));
                sj.add(String.valueOf(breakend.incompleteFragmentCount()));

                sj.add(String.valueOf(breakend.calcQual()));

                String facingBreakendIds = breakend.facingBreakends().stream().map(x -> String.valueOf(x.id())).collect(Collectors.joining(ITEM_DELIM));
                sj.add(facingBreakendIds);

                sj.add(AlternativeAlignment.toVcfTag(breakend.alternativeAlignments()));

                InsertionType insertionType = getInsertionType(breakend, assemblyAlignment.assemblies(), assemblyAlignment.breakends());
                sj.add(insertionType.toString());

                sj.add(String.valueOf(breakend.uniqueFragmentPositionCount()));

                String assemblyMatchStr = getClosestAssembly(breakend, assemblyAlignment.assemblies(), closestAssemblyMap, true);
                sj.add(assemblyMatchStr);

                sj.add(String.valueOf(breakend.nonPrimaryAssemblyFragmentCount()));

                mWriter.write(sj.toString());
                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write breakends: {}", e.toString());
        }
    }

    private static InsertionType getInsertionType(
            final Breakend breakend, final List<JunctionAssembly> assemblies, final List<Breakend> breakends)
    {
        boolean phasedWithLine = assemblies.stream().anyMatch(x -> x.hasLineSequence());

        if(!phasedWithLine)
            return InsertionType.NONE;

        if(isMobileLineElement(breakend.Orient, breakend.InsertedBases))
            return InsertionType.LINE;

        for(Breakend otherBreakend : breakends)
        {
            if(otherBreakend == breakend)
                continue;

            if(withinLineProximity(breakend.Position, otherBreakend.Position, breakend.Orient, otherBreakend.Orient)
            && isMobileLineElement(otherBreakend.Orient, otherBreakend.InsertedBases))
            {
                return InsertionType.NEAR;
            }
        }

        return InsertionType.NONE;
    }

    private static String getClosestAssembly(
            final Breakend breakend, final List<JunctionAssembly> assemblies, final Map<Breakend,String> closestAssemblyMap, boolean checkMate)
    {
        // find an assembly with close to matching coords, otherwise use the coords from a mate breakend, and cache the results
        if(closestAssemblyMap.containsKey(breakend))
            return closestAssemblyMap.get(breakend);

        for(JunctionAssembly assembly : assemblies)
        {
            Junction junction = assembly.junction();

            if(junction.Chromosome.equals(breakend.Chromosome) && junction.Orient == breakend.Orient
            && abs(junction.Position - breakend.Position) < 100)
            {
                String junctionCoords = junction.coordsTyped();

                if(!breakend.isSingle())
                    closestAssemblyMap.put(breakend, junctionCoords);

                return junctionCoords;
            }
        }

        if(breakend.isSingle())
            return "";

        closestAssemblyMap.put(breakend, "");

        if(!checkMate)
            return "";

        Breakend otherBreakend = breakend.otherBreakend();
        String otherBreakendCoords = closestAssemblyMap.get(otherBreakend);

        if(otherBreakendCoords == null)
        {
            otherBreakendCoords = getClosestAssembly(otherBreakend, assemblies, closestAssemblyMap, false);
            closestAssemblyMap.put(otherBreakend, otherBreakendCoords);
        }

        return !otherBreakendCoords.isEmpty() ? otherBreakendCoords + OTHER_BREAKEND_COORDS : "";
    }

    private static final String OTHER_BREAKEND_COORDS = "_OTHER";
}
