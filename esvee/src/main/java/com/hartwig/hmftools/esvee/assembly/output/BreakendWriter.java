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
import com.hartwig.hmftools.esvee.utils.TruthsetAnnotation;

public class BreakendWriter
{
    private final AssemblyConfig mConfig;

    private final BufferedWriter mWriter;
    private final TruthsetAnnotation mTruthsetAnnotation;

    public BreakendWriter(final AssemblyConfig config)
    {
        mConfig = config;
        mTruthsetAnnotation = new TruthsetAnnotation(mConfig.TruthsetFile);

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
            sj.add("AssemblyInfo");

            sj.add("Type").add("Chromosome").add("Position").add("Orientation");

            sj.add("MateChr").add("MatePos").add("MateOrient");
            sj.add("InsertedBases").add("Homology").add("ConfidenceInterval").add("InexactOffset");

            sj.add("SplitFrags").add("DiscFrags");

            sj.add("Filters");

            sj.add("AnchorLength");
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
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            String assemblyInfo = assemblyAlignment.info();

            for(Breakend breakend : assemblyAlignment.breakends())
            {
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

                sj.add(String.valueOf(breakend.anchorLength()));

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
