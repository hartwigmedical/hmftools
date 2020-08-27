package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.FragmentType.ALT;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.isofox.common.FragmentType.UNSPLICED;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.isofox.adjusts.FragmentSizeCalcs;
import com.sun.org.apache.bcel.internal.generic.DUP;

import org.immutables.value.Value;

@Value.Immutable
public abstract class SummaryStats
{
    public abstract int totalFragments();
    public abstract int duplicateFragments();
    public abstract double splicedFragmentPerc();
    public abstract double unsplicedFragmentPerc();
    public abstract double altFragmentPerc();
    public abstract double chimericFragmentPerc();

    public abstract int readLength();

    // 5th, 50th and 95th percentile intronic fragment length
    public abstract double fragmentLength5thPercent();
    public abstract double fragmentLength50thPercent();
    public abstract double fragmentLength95thPercent();

    // proportion of fragments in 7 highly expressed genes
    public abstract double enrichedGenePercent();

    // Median GC (excluding 7 highly expressed genes)
    public abstract double medianGCRatio();

    public static SummaryStats createSummaryStats(
            final int[] totalCounts, int enrichedGeneFragCount,
            double medianGCRatio, final List<int[]> fragmentLengths, int maxReadLength)
    {
        int totalFragments = totalCounts[typeAsInt(TOTAL)];
        int totalDuplicates = totalCounts[typeAsInt(DUPLICATE)];
        double enrichedGenePercent = totalFragments > 0 ? enrichedGeneFragCount / (double)totalFragments : 0;

        final List<Double> fragLengths = FragmentSizeCalcs.calcPercentileData(fragmentLengths, Lists.newArrayList(0.05, 0.5, 0.95));

        return ImmutableSummaryStats.builder()
                .totalFragments(totalFragments)
                .duplicateFragments(totalDuplicates)
                .splicedFragmentPerc(totalCounts[typeAsInt(TRANS_SUPPORTING)] / (double)totalFragments)
                .unsplicedFragmentPerc(totalCounts[typeAsInt(UNSPLICED)] / (double)totalFragments)
                .altFragmentPerc(totalCounts[typeAsInt(ALT)] / (double)totalFragments)
                .chimericFragmentPerc(totalCounts[typeAsInt(CHIMERIC)] / (double)totalFragments)
                .readLength(maxReadLength)
                .fragmentLength5thPercent(!fragLengths.isEmpty() ? fragLengths.get(0) : 0)
                .fragmentLength50thPercent(!fragLengths.isEmpty() ? fragLengths.get(1) : 0)
                .fragmentLength95thPercent(!fragLengths.isEmpty() ? fragLengths.get(2) : 0)
                .enrichedGenePercent(enrichedGenePercent)
                .medianGCRatio(medianGCRatio)
                .build();
    }

    public static String csvHeader()
    {
        return "SampleId,TotalFragments,DuplicateFragments"
                + ",SplicedFragmentPerc,UnsplicedFragmentPerc,AltFragmentPerc,ChimericFragmentPerc"
                + ",ReadLength,FragLength5th,FragLength50th,FragLength95th"
                + ",EnrichedGenePercent,MedianGCRatio";
    }

    public String toCsv(final String sampleId)
    {
        return new StringJoiner(DELIMITER)
                .add(sampleId)
                .add(String.valueOf(totalFragments()))
                .add(String.valueOf(duplicateFragments()))
                .add(String.format("%.3f", splicedFragmentPerc()))
                .add(String.format("%.3f", unsplicedFragmentPerc()))
                .add(String.format("%.3f", altFragmentPerc()))
                .add(String.format("%.3f", chimericFragmentPerc()))
                .add(String.valueOf(readLength()))
                .add(String.format("%.0f", fragmentLength5thPercent()))
                .add(String.format("%.0f", fragmentLength50thPercent()))
                .add(String.format("%.0f", fragmentLength95thPercent()))
                .add(String.format("%.3f", enrichedGenePercent()))
                .add(String.format("%.3f", medianGCRatio()))
                .toString();
    }

    public static SummaryStats fromCsv(final String input)
    {
        final String items[] = input.split(DELIMITER);
        int index = 1;

        return ImmutableSummaryStats.builder()
                .totalFragments(Integer.parseInt(items[index++]))
                .duplicateFragments(Integer.parseInt(items[index++]))
                .splicedFragmentPerc(Double.parseDouble(items[index++]))
                .unsplicedFragmentPerc(Double.parseDouble(items[index++]))
                .altFragmentPerc(Double.parseDouble(items[index++]))
                .chimericFragmentPerc(Double.parseDouble(items[index++]))
                .readLength(Integer.parseInt(items[index++]))
                .fragmentLength5thPercent(Double.parseDouble(items[index++]))
                .fragmentLength50thPercent(Double.parseDouble(items[index++]))
                .fragmentLength95thPercent(Double.parseDouble(items[index++]))
                .enrichedGenePercent(Double.parseDouble(items[index++]))
                .medianGCRatio(Double.parseDouble(items[index++]))
                .build();
    }

    public static SummaryStats loadFile(final Path filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(filename);
            if(lines.size() != 2)
            {
                ISF_LOGGER.error("failed to load file({}: invalid line count({})", filename.toString(), lines.size());
                return null;
            }

            return fromCsv(lines.get(1));
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load summary data file({}): {}", filename.toString(), e.toString());
            return null;
        }
    }

}
