package com.hartwig.hmftools.sage.coverage;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sage.SageCommon.DELIM;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

public final class ExonMedianDepth
{
    public static void write(final String filename, final List<GeneCoverage> geneCoverages)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);
            writer.write(header());
            writer.newLine();

            for(GeneCoverage geneCoverage : geneCoverages)
            {
                for(ExonCoverage exonCoverage : geneCoverage.exonCoverage())
                {
                    StringJoiner data = new StringJoiner(DELIM);
                    data.add(geneCoverage.geneName());
                    data.add(exonCoverage.chromosome());
                    data.add(String.valueOf(exonCoverage.start()));
                    data.add(String.valueOf(exonCoverage.end()));
                    data.add(String.valueOf(exonCoverage.exonRank()));
                    data.add(format("%.0f", exonCoverage.medianDepth()));
                    writer.write(data.toString());
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write exon median coverage: {}", e.toString());
        }
    }

    private static String header()
    {
        StringJoiner header = new StringJoiner(DELIM);
        header.add("gene");
        header.add("chromosome");
        header.add("posStart");
        header.add("posEnd");
        header.add("exonRank");
        header.add("medianDepth");
        return header.toString();
    }
}
