package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.patientreporter.slicing.GenomeRegion;

import org.jetbrains.annotations.NotNull;

final class CopyNumberAnalyser {

    private static final int NORMAL_COPYNUMBER = 2;

    private CopyNumberAnalyser() {
    }

    @NotNull
    static Map<GenomeRegion, CopyNumberStats> run(@NotNull Collection<GenomeRegion> regions,
            @NotNull Collection<CopyNumber> copyNumbers) {
        final Map<GenomeRegion, CopyNumberStats> stats = Maps.newHashMap();

        final Multimap<String, CopyNumber> copyNumberMap = toMultiMap(copyNumbers);
        for (final GenomeRegion region : regions) {
            stats.put(region, analyzeRegion(region, Lists.newArrayList(copyNumberMap.get(region.chromosome()))));
        }
        return stats;
    }

    @NotNull
    private static Multimap<String, CopyNumber> toMultiMap(@NotNull final Collection<CopyNumber> copyNumbers) {
        final Multimap<String, CopyNumber> copyNumberMultimap = ArrayListMultimap.create();
        for (final CopyNumber copyNumber : copyNumbers) {
            copyNumberMultimap.put(copyNumber.chromosome(), copyNumber);
        }
        return copyNumberMultimap;
    }

    @NotNull
    private static CopyNumberStats analyzeRegion(@NotNull final GenomeRegion region,
            @NotNull final List<CopyNumber> copyNumbers) {
        int minCopyNumber = Integer.MAX_VALUE;
        int maxCopyNumber = Integer.MIN_VALUE;
        long totalCopyNumber = NORMAL_COPYNUMBER * region.bases();

        long totalOverlap = 0L;
        for (final CopyNumber copyNumber : copyNumbers) {
            final long overlap = overlappingBases(region, copyNumber);
            if (overlap > 0) {
                minCopyNumber = Math.min(copyNumber.value(), minCopyNumber);
                maxCopyNumber = Math.max(copyNumber.value(), maxCopyNumber);
                totalCopyNumber += (overlap * (copyNumber.value() - NORMAL_COPYNUMBER));
            }
            totalOverlap += overlap;
        }

        if (totalOverlap == 0) {
            minCopyNumber = NORMAL_COPYNUMBER;
            maxCopyNumber = NORMAL_COPYNUMBER;
        } else if (totalOverlap < region.bases()) {
            minCopyNumber = Math.min(NORMAL_COPYNUMBER, minCopyNumber);
            maxCopyNumber = Math.max(NORMAL_COPYNUMBER, maxCopyNumber);
        }

        return new CopyNumberStats(minCopyNumber, maxCopyNumber, ((double) totalCopyNumber) / region.bases());
    }

    private static long overlappingBases(@NotNull final GenomeRegion region, @NotNull final CopyNumber copyNumber) {
        long minEnd = Math.min(copyNumber.end(), region.end());
        long maxStart = Math.max(copyNumber.start(), region.start());
        return Math.max(0, 1 + minEnd - maxStart);
    }
}
