package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.ChromosomeArm;

public record ChromosomeArmCopyNumber(HumanChromosome chromosome,
                                      ChromosomeArm arm,
                                      double meanCopyNumber,
                                      double medianCopyNumber,
                                      double minCopyNumber,
                                      double maxCopyNumber)
{
    static final DecimalFormat FORMAT = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.ENGLISH));

    public static String tsvFileHeader()
    {
        return new StringJoiner("\t")
                .add("chromosome")
                .add("arm")
                .add("meanCopyNumber")
                .add("medianCopyNumber")
                .add("minCopyNumber")
                .add("maxCopyNumber")
                .toString();
    }

    public String toTSV()
    {
        return new StringJoiner(TSV_DELIM)
                .add(chromosome.toString())
                .add(ChromosomeArm.asStr(arm))
                .add(FORMAT.format(meanCopyNumber))
                .add(FORMAT.format(medianCopyNumber))
                .add(FORMAT.format(minCopyNumber))
                .add(FORMAT.format(maxCopyNumber))
                .toString();
    }

    public boolean includeInReport()
    {
        if(chromosome.hasShortArm())
        {
            return arm == ChromosomeArm.Q_ARM;
        }
        return true;
    }
}
