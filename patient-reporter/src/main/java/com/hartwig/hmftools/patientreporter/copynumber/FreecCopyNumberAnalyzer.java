package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.slicing.HmfSlicer;

import org.jetbrains.annotations.NotNull;

public final class FreecCopyNumberAnalyzer {

    private static final int NORMAL_COPYNUMBER = 2;
    private static final int MIN_CNV_FOR_GAIN = 4;
    private static final int MAX_CNV_FOR_LOSS = 0;

    @NotNull
    private final Set<HmfGenomeRegion> regions;

    public static FreecCopyNumberAnalyzer fromHmfSlicingRegion(@NotNull final HmfSlicer hmfSlicingRegion) {
        return new FreecCopyNumberAnalyzer(Sets.newHashSet(hmfSlicingRegion.hmfRegions()));
    }

    private FreecCopyNumberAnalyzer(@NotNull final Set<HmfGenomeRegion> regions) {
        this.regions = regions;
    }

    @NotNull
    public CopyNumberAnalysis run(@NotNull final List<CopyNumber> copyNumbers) {
        final Map<HmfGenomeRegion, CopyNumberStats> stats = Maps.newHashMap();

        final Multimap<String, CopyNumber> copyNumberPerChromosome = toChromosomeMultiMap(copyNumbers);
        for (final HmfGenomeRegion region : regions) {
            stats.put(region,
                    analyzeRegion(region, Lists.newArrayList(copyNumberPerChromosome.get(region.chromosome()))));
        }
        final List<CopyNumberReport> reports = Lists.newArrayList();
        for (final Map.Entry<HmfGenomeRegion, CopyNumberStats> stat : stats.entrySet()) {
            final int relevantCNV = stat.getValue().min();
            if (relevantCNV >= MIN_CNV_FOR_GAIN || relevantCNV <= MAX_CNV_FOR_LOSS) {
                final HmfGenomeRegion region = stat.getKey();
                reports.add(ImmutableCopyNumberReport.builder()
                        .chromosome(stat.getKey().chromosome())
                        .chromosomeBand(region.chromosomeBand())
                        .gene(region.gene())
                        .copyNumber(relevantCNV)
                        .type(CopyNumberReportType.resolveType(relevantCNV))
                        .build());
            }
        }
        Collections.sort(reports);
        return new CopyNumberAnalysis(stats.size(), reports);
    }

    @NotNull
    private static Multimap<String, CopyNumber> toChromosomeMultiMap(
            @NotNull final Collection<CopyNumber> copyNumbers) {
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
            final long overlap = region.overlappingBases(copyNumber);
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
}
