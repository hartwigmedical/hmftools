package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.patientreporter.slicing.GenomeRegion;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CopyNumberAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberAnalyzer.class);

    private static final int NORMAL_COPYNUMBER = 2;
    private static final int MIN_CNV_FOR_GAIN = 4;
    private static final int MAX_CNV_FOR_LOSS = 0;

    @NotNull
    private final Map<GenomeRegion, HMFSlicingAnnotation> annotations;

    public static CopyNumberAnalyzer fromHmfSlicingRegion(@NotNull final Slicer hmfSlicingRegion) {
        final Map<GenomeRegion, HMFSlicingAnnotation> annotations = Maps.newHashMap();
        for (final GenomeRegion region : hmfSlicingRegion.regions()) {
            final HMFSlicingAnnotation annotation = HMFSlicingAnnotation.fromGenomeRegion(region);
            if (annotation == null) {
                LOGGER.warn("Could not extract annotation from hmf slicing region: " + region);
            } else {
                annotations.put(region, annotation);
            }
        }
        return new CopyNumberAnalyzer(annotations);
    }

    @VisibleForTesting
    CopyNumberAnalyzer(@NotNull final Map<GenomeRegion, HMFSlicingAnnotation> annotations) {
        this.annotations = annotations;
    }

    @NotNull
    public CopyNumberAnalysis run(@NotNull List<CopyNumber> copyNumbers) {
        final Map<GenomeRegion, CopyNumberStats> stats = Maps.newHashMap();

        final Multimap<String, CopyNumber> copyNumberPerChromosome = toChromosomeMultiMap(copyNumbers);
        for (final GenomeRegion region : annotations.keySet()) {
            stats.put(region,
                    analyzeRegion(region, Lists.newArrayList(copyNumberPerChromosome.get(region.chromosome()))));
        }
        final List<CopyNumberReport> reports = Lists.newArrayList();
        for (final Map.Entry<GenomeRegion, CopyNumberStats> stat : stats.entrySet()) {
            final int relevantCNV = stat.getValue().min();

            if (relevantCNV >= MIN_CNV_FOR_GAIN || relevantCNV <= MAX_CNV_FOR_LOSS) {
                final HMFSlicingAnnotation annotation = annotations.get(stat.getKey());
                reports.add(new CopyNumberReport.Builder().gene(annotation.gene()).
                        transcript(annotation.transcript()).copyNumber(relevantCNV).build());
            }
        }
        return new CopyNumberAnalysis(stats, reports);
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
