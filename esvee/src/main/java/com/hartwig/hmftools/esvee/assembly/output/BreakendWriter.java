package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.filters.FilterType.filtersAsStr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.AssemblyConfig;
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

            // sj.add("Id");
            sj.add("AssemblyId");
            sj.add("AssemblyInfo");

            sj.add("Type").add("Chromosome").add("Position").add("Orientation");

            sj.add("MateChr").add("MatePos").add("MateOrient");
            sj.add("InsertedBases").add("Homology").add("ConfidenceInterval").add("InexactOffset");

            sj.add("SplitFrags").add("DiscFrags");

            sj.add("Filters");

            sj.add("SequenceLength");
            sj.add("SegmentCount");
            sj.add("SequenceIndex");
            sj.add("SegmentOrientation");
            sj.add("SegmentIndex");
            sj.add("AlignedBases");
            sj.add("MapQual");
            sj.add("Score");
            sj.add("RepeatTrimLength");
            sj.add("SplitFragments");
            sj.add("RefSplitFragments");
            sj.add("DiscFragments");
            sj.add("RefDiscFragments");
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
                sj.add(String.valueOf(assemblyAlignment.id()));
                sj.add(assemblyInfo);
                sj.add(String.valueOf(breakend.svType()));
                sj.add(breakend.Chromosome);
                sj.add(String.valueOf(breakend.Position));
                sj.add(String.valueOf(breakend.Orientation));

                if(breakend.otherBreakend() != null)
                {
                    sj.add(breakend.otherBreakend().Chromosome);
                    sj.add(String.valueOf(breakend.otherBreakend().Position));
                    sj.add(String.valueOf(breakend.otherBreakend().Orientation));

                }
                else
                {
                    sj.add("").add("").add("");
                }

                sj.add(breakend.InsertedBases);

                if(breakend.Homology != null)
                {
                    sj.add(breakend.Homology.Homology);
                    sj.add(format("%d-%d", breakend.Homology.ExactStart, breakend.Homology.ExactEnd));
                    sj.add(format("%d-%d", breakend.Homology.InexactStart, breakend.Homology.InexactEnd));
                }
                else
                {
                    sj.add("").add("0-0").add("0-0");
                }

                int totalSplitFrags = breakend.sampleSupport().stream().mapToInt(x -> x.SplitFragments).sum();
                int totalDiscFrags = breakend.sampleSupport().stream().mapToInt(x -> x.DiscordantFragments).sum();

                sj.add(String.valueOf(totalSplitFrags));
                sj.add(String.valueOf(totalDiscFrags));

                sj.add(filtersAsStr(breakend.filters()));

                // sj.add(String.valueOf(breakend.anchorLength()));
                sj.add(String.valueOf(assemblyAlignment.fullSequenceLength()));

                // for now just the first segment - no showing branching or duplicates
                BreakendSegment segment = breakend.segments().get(0);
                sj.add(String.valueOf(breakend.segments().size()));
                sj.add(String.valueOf(segment.SequenceIndex));
                sj.add(String.valueOf(segment.Orientation));
                sj.add(String.valueOf(segment.Index));
                sj.add(String.valueOf(segment.Alignment.alignedBases()));
                sj.add(String.valueOf(segment.Alignment.MapQual));
                sj.add(String.valueOf(segment.Alignment.Score));
                sj.add(String.valueOf(segment.Alignment.repeatTrimmedLength()));

                int tumorCount = mConfig.TumorIds.size();

                int splitFrags = 0;
                int refSplitFrags = 0;
                int discFrags = 0;
                int refDiscFrags = 0;

                for(int i = 0; i < breakend.sampleSupport().size(); ++i)
                {
                    splitFrags += breakend.sampleSupport().get(i).SplitFragments;
                    discFrags += breakend.sampleSupport().get(i).DiscordantFragments;

                    if(i >= tumorCount)
                    {
                        refSplitFrags += breakend.sampleSupport().get(i).SplitFragments;
                        refDiscFrags += breakend.sampleSupport().get(i).DiscordantFragments;
                    }
                }

                sj.add(String.valueOf(splitFrags));
                sj.add(String.valueOf(refSplitFrags));
                sj.add(String.valueOf(discFrags));
                sj.add(String.valueOf(refDiscFrags));

                String altAlignmentsStr = breakend.alternativeAlignments().stream()
                        .map(x -> x.altAlignmentStr()).collect(Collectors.joining(";"));
                sj.add(altAlignmentsStr);

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
