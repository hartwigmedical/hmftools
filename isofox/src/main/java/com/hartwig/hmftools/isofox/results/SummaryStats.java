package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.common.rna.RnaQcFilter.qcFiltersToString;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.FragmentType.ALT;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.FORWARD_STRAND;
import static com.hartwig.hmftools.isofox.common.FragmentType.REVERSE_STRAND;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.isofox.common.FragmentType.UNSPLICED;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.ImmutableRnaStatistics;
import com.hartwig.hmftools.common.rna.RnaQcFilter;
import com.hartwig.hmftools.common.rna.RnaStatisticFile;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.adjusts.FragmentSizeCalcs;
import com.hartwig.hmftools.isofox.common.FragmentTypeCounts;

public class SummaryStats
{
    public static RnaStatistics createSummaryStats(
            final FragmentTypeCounts fragmentTypeCounts, long enrichedGeneFragCount, int spliceGeneCount,
            double medianGCRatio, final List<FragmentSize> fragmentLengths, int maxReadLength,
            int lowCoverageThreshold, int splicedGeneThreshold)
    {
        long totalFragments = fragmentTypeCounts.typeCount(TOTAL);
        long duplicateFragments = fragmentTypeCounts.typeCount(DUPLICATE);

        double totalFragmentsDenom = totalFragments;

        double enrichedGenePercent = totalFragments > 0 ? enrichedGeneFragCount / totalFragmentsDenom : 0;

        long fowardFrags = fragmentTypeCounts.typeCount(FORWARD_STRAND);
        double totalStrandFrags = fowardFrags + fragmentTypeCounts.typeCount(REVERSE_STRAND);
        double forwardStrandPerc = totalStrandFrags > 0 ? fowardFrags / totalStrandFrags : 0;

        final List<Double> fragLengths = FragmentSizeCalcs.calcPercentileData(fragmentLengths, Lists.newArrayList(0.05, 0.5, 0.95));

        List<RnaQcFilter> qcFilters = RnaStatisticFile.calcQcStatus(
                totalFragments, duplicateFragments, spliceGeneCount, lowCoverageThreshold, splicedGeneThreshold);

        return ImmutableRnaStatistics.builder()
                .qcStatus(qcFilters)
                .totalFragments(totalFragments)
                .duplicateFragments(duplicateFragments)
                .splicedFragmentPerc(fragmentTypeCounts.typeCount(TRANS_SUPPORTING) / totalFragmentsDenom)
                .unsplicedFragmentPerc(fragmentTypeCounts.typeCount(UNSPLICED) / totalFragmentsDenom)
                .altFragmentPerc(fragmentTypeCounts.typeCount(ALT) / totalFragmentsDenom)
                .chimericFragmentPerc(fragmentTypeCounts.typeCount(CHIMERIC) / totalFragmentsDenom)
                .splicedGeneCount(spliceGeneCount)
                .readLength(maxReadLength)
                .fragmentLength5thPercent(!fragLengths.isEmpty() ? fragLengths.get(0) : 0)
                .fragmentLength50thPercent(!fragLengths.isEmpty() ? fragLengths.get(1) : 0)
                .fragmentLength95thPercent(!fragLengths.isEmpty() ? fragLengths.get(2) : 0)
                .enrichedGenePercent(enrichedGenePercent)
                .medianGCRatio(medianGCRatio)
                .forwardStrandPercent(forwardStrandPerc)
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

            return RnaStatisticFile.fromLines(lines);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load summary data file({}): {}", filename.toString(), e.toString());
            return null;
        }
    }

}
