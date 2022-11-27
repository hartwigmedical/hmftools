package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.FragmentType.ALT;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.isofox.common.FragmentType.UNSPLICED;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.ImmutableRnaStatistics;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.adjusts.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.common.FragmentTypeCounts;

public class SummaryStats
{
    public static RnaStatistics createSummaryStats(
            final FragmentTypeCounts fragmentTypeCounts, long enrichedGeneFragCount,
            double medianGCRatio, final List<FragmentSize> fragmentLengths, int maxReadLength)
    {
        double totalFragments = fragmentTypeCounts.typeCount(TOTAL);
        double enrichedGenePercent = totalFragments > 0 ? enrichedGeneFragCount / (double)totalFragments : 0;

        final List<Double> fragLengths = FragmentSizeCalcs.calcPercentileData(fragmentLengths, Lists.newArrayList(0.05, 0.5, 0.95));

        return ImmutableRnaStatistics.builder()
                .totalFragments(fragmentTypeCounts.typeCount(TOTAL))
                .duplicateFragments(fragmentTypeCounts.typeCount(DUPLICATE))
                .splicedFragmentPerc(fragmentTypeCounts.typeCount(TRANS_SUPPORTING) / totalFragments)
                .unsplicedFragmentPerc(fragmentTypeCounts.typeCount(UNSPLICED) / totalFragments)
                .altFragmentPerc(fragmentTypeCounts.typeCount(ALT) / totalFragments)
                .chimericFragmentPerc(fragmentTypeCounts.typeCount(CHIMERIC) / totalFragments)
                .readLength(maxReadLength)
                .fragmentLength5thPercent(!fragLengths.isEmpty() ? fragLengths.get(0) : 0)
                .fragmentLength50thPercent(!fragLengths.isEmpty() ? fragLengths.get(1) : 0)
                .fragmentLength95thPercent(!fragLengths.isEmpty() ? fragLengths.get(2) : 0)
                .enrichedGenePercent(enrichedGenePercent)
                .medianGCRatio(medianGCRatio)
                .build();
    }

    public static RnaStatistics loadFile(final Path filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(filename);
            if(lines.size() != 2)
            {
                ISF_LOGGER.error("failed to load file({}: invalid line count({})", filename.toString(), lines.size());
                return null;
            }

            return RnaStatistics.fromCsv(lines.get(1));
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load summary data file({}): {}", filename.toString(), e.toString());
            return null;
        }
    }

}
