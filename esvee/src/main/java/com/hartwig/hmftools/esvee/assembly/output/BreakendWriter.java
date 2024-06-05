package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FilterType.filtersAsStr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.alignment.AlternativeAlignment;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.alignment.Breakend;
import com.hartwig.hmftools.esvee.alignment.BreakendSegment;
import com.hartwig.hmftools.esvee.utils.TruthsetAnnotation;

public class BreakendWriter
{
    private final AssemblyConfig mConfig;

    private final BufferedWriter mWriter;
    private final TruthsetAnnotation mTruthsetAnnotation;

    public BreakendWriter(final AssemblyConfig config, final TruthsetAnnotation truthsetAnnotation)
    {
        mConfig = config;
        mTruthsetAnnotation = truthsetAnnotation;

        mWriter = initialiseWriter();
    }

    public void close() { closeBufferedWriter(mWriter);}

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
            sj.add("AssemblyId");
            sj.add("AssemblyInfo");

            sj.add("Type").add("Chromosome").add("Position").add("Orientation");

            sj.add("MateChr").add("MatePos").add("MateOrient").add("Length");
            sj.add("InsertedBases").add("Homology").add("ConfidenceInterval").add("InexactOffset");

            sj.add("Filters");

            sj.add("Qual");
            sj.add("SplitFragments");
            sj.add("RefSplitFragments");
            sj.add("DiscFragments");
            sj.add("RefDiscFragments");
            sj.add("ForwardReads");
            sj.add("ReverseReads");

            sj.add("SequenceLength");
            sj.add("SegmentCount");
            sj.add("SequenceIndex");
            sj.add("SegmentOrientation");
            sj.add("SegmentIndex");
            sj.add("AlignedBases");
            sj.add("MapQual");
            sj.add("Score");
            sj.add("AdjAlignedBases");
            sj.add("AvgFragmentLength");
            sj.add("BreakendQual");

            sj.add("FacingBreakendIds");

            sj.add("AltAlignments");

            if(mTruthsetAnnotation.enabled())
                sj.add(TruthsetAnnotation.tsvHeader());

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

            for(Breakend breakend : assemblyAlignment.breakends())
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);

                sj.add(String.valueOf(breakend.id()));
                sj.add(String.valueOf(assemblyAlignment.assemblies().get(0).phaseGroup().id()));
                sj.add(String.valueOf(assemblyAlignment.id()));
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

                if(breakend.Homology != null)
                {
                    sj.add(breakend.Homology.Homology);
                    sj.add(format("%d,%d", breakend.Homology.ExactStart, breakend.Homology.ExactEnd));
                    sj.add(format("%d,%d", breakend.Homology.InexactStart, breakend.Homology.InexactEnd));
                }
                else
                {
                    sj.add("").add("0,0").add("0,0");
                }

                sj.add(filtersAsStr(breakend.filters()));

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
                sj.add(String.valueOf(segment.SequenceIndex));
                sj.add(String.valueOf(segment.Orient));
                sj.add(String.valueOf(segment.Index));
                sj.add(String.valueOf(segment.Alignment.alignedBases()));
                sj.add(String.valueOf(segment.Alignment.MapQual));
                sj.add(String.valueOf(segment.Alignment.Score));
                sj.add(String.valueOf(segment.Alignment.adjustedAlignment()));
                sj.add(String.valueOf(breakend.averageFragmentLength()));

                sj.add(String.valueOf(breakend.calcQual()));

                String facingBreakendIds = breakend.facingBreakends().stream().map(x -> String.valueOf(x.id())).collect(Collectors.joining(ITEM_DELIM));
                sj.add(facingBreakendIds);

                sj.add(AlternativeAlignment.toVcfTag(breakend.alternativeAlignments()));

                if(mTruthsetAnnotation.enabled())
                    sj.add(mTruthsetAnnotation.findTruthsetAnnotation(breakend));

                mWriter.write(sj.toString());
                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write breakends: {}", e.toString());
        }
    }
}
