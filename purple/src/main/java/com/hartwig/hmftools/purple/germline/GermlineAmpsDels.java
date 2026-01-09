package com.hartwig.hmftools.purple.germline;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.ANY;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.VARIANT_NOT_LOST;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.WILDTYPE_LOST;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.ReportedStatus.NONE;
import static com.hartwig.hmftools.common.purple.ReportedStatus.REPORTED;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_CN_CONSISTENCY_MACN_PERC;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_CN_CONSISTENCY_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_COHORT_FREQ;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_GENE_BUFFER;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_NORMAL_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_REGION_MATCH_BUFFER;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_AMP_DEL_REGION_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.drivers.DeletionDrivers.MAX_COPY_NUMBER_DEL;
import static com.hartwig.hmftools.purple.germline.GermlineDeletionUtils.findMatchingCopyNumber;
import static com.hartwig.hmftools.purple.germline.GermlineDeletionUtils.findOverlappingExons;
import static com.hartwig.hmftools.purple.germline.GermlineDeletionUtils.findOverlappingGenes;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineDetectionMethod;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.jetbrains.annotations.Nullable;

public class GermlineAmpsDels
{
    static final String FILTER_CN_INCONSISTENCY = "INCONSISTENT_CN";
    static final String FILTER_COHORT_FREQ = "COHORT_FREQ";
    static final String FILTER_REGION_LENGTH = "MIN_LENGTH";

    public interface GeneDataSupplier
    {
        List<GeneData> getGeneData(String chromosome);

        TranscriptData getTranscriptData(String geneId);
    }

    private record EnsemblDataSupplier(EnsemblDataCache mCache) implements GeneDataSupplier
    {
        @Override
        public List<GeneData> getGeneData(final String chromosome)
        {
            return mCache.getChrGeneDataMap().getOrDefault(chromosome, List.of());
        }

        @Override
        public TranscriptData getTranscriptData(final String geneId)
        {
            return mCache.getTranscriptData(geneId, "");
        }
    }

    private final Map<String, DriverGene> mDriverGenes;
    private final GeneDataSupplier mGeneDataCache;
    private final GermlineDeletionFrequency mCohortFrequency;

    private final List<GermlineAmpDel> mEvents;
    private final List<DriverCatalog> mDrivers;

    public GermlineAmpsDels(
            final Map<String, DriverGene> driverGenes, final EnsemblDataCache geneDataCache,
            final GermlineDeletionFrequency cohortFrequency)
    {
        this(driverGenes, new EnsemblDataSupplier(geneDataCache), cohortFrequency);
    }

    public GermlineAmpsDels(final Map<String, DriverGene> driverGenes, final GeneDataSupplier geneDataCache,
            final GermlineDeletionFrequency cohortFrequency)
    {
        mGeneDataCache = geneDataCache;
        mDriverGenes = driverGenes;
        mCohortFrequency = cohortFrequency;

        mEvents = Lists.newArrayList();
        mDrivers = Lists.newArrayList();
    }

    public List<GermlineAmpDel> getEvents()
    {
        return mEvents;
    }

    public List<DriverCatalog> getDrivers()
    {
        return mDrivers;
    }

    public void findEvents(
            final List<PurpleCopyNumber> copyNumbers, final List<ObservedRegion> fittedRegions, final List<StructuralVariant> germlineSVs)
    {
        for(int i = 0; i < fittedRegions.size(); ++i)
        {
            ObservedRegion region = fittedRegions.get(i);
            if(region.germlineStatus() != HOM_DELETION && region.germlineStatus() != HET_DELETION)
            {
                continue;
            }
            PurpleCopyNumber overlappingCopyNumberWithLowestMajorAlleleCopyNumber = findMatchingCopyNumber(region, copyNumbers);
            ObservedRegion nextRegion = i < fittedRegions.size() - 1 ? fittedRegions.get(i + 1) : null;
            if(nextRegion != null && !nextRegion.chromosome().equals(region.chromosome()))
            {
                nextRegion = null;
            }
            MatchedStructuralVariant[] matchingSVs = findMatchingGermlineSVs(region, nextRegion, germlineSVs);
            findOverlappingDriverGene(region, overlappingCopyNumberWithLowestMajorAlleleCopyNumber, matchingSVs);
        }
    }

    private static class MatchedStructuralVariant
    {
        public final StructuralVariant Variant;
        public final boolean IsStart;

        public MatchedStructuralVariant(final StructuralVariant variant, final boolean isStart)
        {
            Variant = variant;
            IsStart = isStart;
        }

        public Integer position()
        {
            return Variant.position(IsStart);
        }
    }

    private MatchedStructuralVariant[] findMatchingGermlineSVs(
            final ObservedRegion region, @Nullable final ObservedRegion nextRegion, final List<StructuralVariant> germlineSVs)
    {
        MatchedStructuralVariant[] matchingSVs = new MatchedStructuralVariant[] { null, null };

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int regionStart;
            int regionEnd;
            byte requiredOrientation = (se == SE_START) ? ORIENT_FWD : ORIENT_REV;

            if(se == SE_START)
            {
                regionStart = region.minStart();
                regionEnd = region.maxStart();
            }
            else
            {
                regionEnd = nextRegion != null ? nextRegion.maxStart() - 1 : region.end();
                regionStart = region.end();
            }

            regionStart -= WINDOW_SIZE;
            regionEnd += WINDOW_SIZE;

            // first check if the matched SV for the start also matches at the end
            if(matchingSVs[SE_START] != null && matchingSVs[SE_START].Variant.end() != null)
            {
                StructuralVariant sv = matchingSVs[SE_START].Variant;

                if(legMatchesRegion(sv.end(), region.chromosome(), regionStart, regionEnd, requiredOrientation))
                {
                    matchingSVs[SE_END] = new MatchedStructuralVariant(sv, false);
                    break;
                }
            }

            for(StructuralVariant variant : germlineSVs)
            {
                for(int varSe = SE_START; varSe <= SE_END; ++varSe)
                {
                    boolean isStart = varSe == SE_START;
                    StructuralVariantLeg leg = isStart ? variant.start() : variant.end();

                    if(leg == null)
                    {
                        continue;
                    }

                    if(!legMatchesRegion(leg, region.chromosome(), regionStart, regionEnd, requiredOrientation))
                    {
                        continue;
                    }

                    if(matchingSVs[se] == null || selectNewBySvType(variant.type(), matchingSVs[se].Variant.type()))
                    {
                        matchingSVs[se] = new MatchedStructuralVariant(variant, isStart);
                    }
                }
            }
        }
        return matchingSVs;
    }

    private static boolean legMatchesRegion(
            final StructuralVariantLeg leg, final String regionChr, int regionStart, int regionEnd, byte requiredOrientation)
    {
        return leg.chromosome().equals(regionChr) && positionWithin(leg.position(), regionStart, regionEnd)
                && leg.orientation() == requiredOrientation;
    }

    private static boolean selectNewBySvType(final StructuralVariantType newType, final StructuralVariantType existingType)
    {
        if(existingType == DEL)
        {
            return false;
        }

        return newType == DEL;
    }

    private static List<String> checkFilters(final ObservedRegion region, final PurpleCopyNumber copyNumber, int cohortFrequency)
    {
        final List<String> filters = Lists.newArrayList();

        if(region.end() - region.start() <= GERMLINE_AMP_DEL_REGION_MIN)
        {
            filters.add(FILTER_REGION_LENGTH);
        }

        if(copyNumber != null)
        {
            double cnInconsistency = region.refNormalisedCopyNumber() - copyNumber.majorAlleleCopyNumber();
            double cnLimit =
                    max(GERMLINE_AMP_DEL_CN_CONSISTENCY_MIN, copyNumber.majorAlleleCopyNumber() * GERMLINE_AMP_DEL_CN_CONSISTENCY_MACN_PERC);

            if(cnInconsistency > cnLimit)
            {
                filters.add(FILTER_CN_INCONSISTENCY);
            }
        }

        if(!filters.contains(FILTER_CN_INCONSISTENCY) && region.observedNormalRatio() > GERMLINE_AMP_DEL_NORMAL_RATIO)
        {
            filters.add(FILTER_CN_INCONSISTENCY);
        }

        if(cohortFrequency >= GERMLINE_AMP_DEL_COHORT_FREQ)
        {
            filters.add(FILTER_COHORT_FREQ);
        }

        return filters;
    }

    private void findOverlappingDriverGene(
            final ObservedRegion region, final PurpleCopyNumber matchedCopyNumber, final MatchedStructuralVariant[] matchingSVs)
    {
        // now find genes
        List<GeneData> geneDataList = mGeneDataCache.getGeneData(region.chromosome());
        if(geneDataList == null)
        {
            return;
        }

        int adjustPosStart = matchingSVs[SE_START] != null ? matchingSVs[SE_START].position() : region.start();
        int adjustPosEnd = matchingSVs[SE_END] != null ? matchingSVs[SE_END].position() : region.end();

        double tumorCopyNumber = region.refNormalisedCopyNumber();
        GermlineStatus tumorStatus = tumorCopyNumber < 0.5 ? HOM_DELETION : HET_DELETION;

        // find overlapping driver genes
        List<GeneData> overlappingGenes =
                findOverlappingGenes(region.chromosome(), adjustPosStart, adjustPosEnd, GERMLINE_AMP_DEL_GENE_BUFFER, geneDataList);

        List<DriverGene> driverGenes = Lists.newArrayList();
        List<TranscriptData> transcripts = Lists.newArrayList();

        for(GeneData geneData : overlappingGenes)
        {
            DriverGene driverGene = mDriverGenes.get(geneData.GeneName);
            if(driverGene != null)
            {
                driverGenes.add(driverGene);
            }
        }

        if(overlappingGenes.isEmpty())
        {
            return;
        }

        int cohortFrequency =
                mCohortFrequency.getRegionFrequency(region.chromosome(), region.start(), region.end(), GERMLINE_AMP_DEL_REGION_MATCH_BUFFER);

        List<String> filters = checkFilters(region, matchedCopyNumber, cohortFrequency);

        List<GeneData> deletedGenes = Lists.newArrayList();
        List<int[]> deletedExonRanges = Lists.newArrayList();

        for(GeneData geneData : overlappingGenes)
        {
            TranscriptData transData = mGeneDataCache.getTranscriptData(geneData.GeneId);
            if(transData == null)
            {
                continue;
            }

            List<ExonData> overlappedExons = findOverlappingExons(transData, adjustPosStart, adjustPosEnd, GERMLINE_AMP_DEL_GENE_BUFFER);
            if(overlappedExons.isEmpty())
            {
                continue;
            }

            PPL_LOGGER.trace("region({}: {}-{}) overlaps gene({}) exons({})",
                    region.chromosome(), region.start(), region.end(), geneData.GeneName, overlappedExons.size());

            int exonRankMin = overlappedExons.stream().mapToInt(x -> x.Rank).min().orElse(0);
            int exonRankMax = overlappedExons.stream().mapToInt(x -> x.Rank).max().orElse(0);

            transcripts.add(transData);
            deletedGenes.add(geneData);
            deletedExonRanges.add(new int[] { exonRankMin, exonRankMax });
        }

        if(transcripts.isEmpty())
        {
            return;
        }

        double germlineCopyNumber = region.observedNormalRatio() * 2;
        String filter = filters.isEmpty() ? PASS_FILTER : String.join(";", filters);
        for(int i = 0; i < deletedGenes.size(); ++i)
        {
            GeneData geneData = deletedGenes.get(i);
            int[] deletedExonRange = deletedExonRanges.get(i);

            ReportedStatus reportedStatus = NONE;

            if(filters.isEmpty())
            {
                for(DriverGene driverGene : driverGenes)
                {
                    if(driverGene.gene().equals(geneData.GeneName) && reportGermlineDeletion(driverGene, tumorStatus))
                    {
                        reportedStatus = REPORTED;
                    }
                }
            }

            mEvents.add(new GermlineAmpDel(
                    geneData.GeneName, region.chromosome(), geneData.KaryotypeBand, adjustPosStart, adjustPosEnd,
                    region.depthWindowCount(), deletedExonRange[0], deletedExonRange[1],
                    GermlineDetectionMethod.SEGMENT, region.germlineStatus(), tumorStatus, germlineCopyNumber, tumorCopyNumber,
                    filter, cohortFrequency, reportedStatus));
        }

        for(TranscriptData transData : transcripts)
        {
            GeneData geneData = overlappingGenes.stream().filter(x -> x.GeneId.equals(transData.GeneId)).findFirst().orElse(null);
            DriverGene driverGene = driverGenes.stream().filter(x -> x.gene().equals(geneData.GeneName)).findFirst().orElse(null);

            if(driverGene == null)
            {
                continue;
            }

            ReportedStatus reportedStatus = NONE;

            if(filters.isEmpty() && reportGermlineDeletion(driverGene, tumorStatus))
            {
                reportedStatus = REPORTED;
            }

            // only create one record even if multiple sections of the gene are deleted
            if(mDrivers.stream().anyMatch(x -> x.gene().equals(geneData.GeneName) && x.transcript().equals(transData.TransName)))
            {
                continue;
            }

            // create a driver record for reportable genes
            DriverCatalog driverCatalog = ImmutableDriverCatalog.builder()
                    .chromosome(region.chromosome())
                    .chromosomeBand(geneData.KaryotypeBand)
                    .gene(geneData.GeneName)
                    .transcript(transData.TransName)
                    .isCanonical(transData.IsCanonical)
                    .driver(DriverType.GERMLINE_DELETION)
                    .category(driverGene.likelihoodType())
                    .driverLikelihood(1)
                    .reportedStatus(reportedStatus)
                    .missense(0)
                    .nonsense(0)
                    .splice(0)
                    .inframe(0)
                    .frameshift(0)
                    .biallelic(region.refNormalisedCopyNumber() < MAX_COPY_NUMBER_DEL)
                    .minCopyNumber(region.refNormalisedCopyNumber())
                    .maxCopyNumber(region.refNormalisedCopyNumber())
                    .likelihoodMethod(LikelihoodMethod.GERMLINE)
                    .build();

            mDrivers.add(driverCatalog);
        }
    }

    private static boolean reportGermlineDeletion(final DriverGene driverGene, final GermlineStatus germlineStatus)
    {
        // check requirements on the germline disruption field: WILDTYPE_LOST - requires an LOH for the deletion to be reportable
        return driverGene.reportGermlineDeletion() == ANY
                || driverGene.reportGermlineDeletion() == VARIANT_NOT_LOST
                || (driverGene.reportGermlineDeletion() == WILDTYPE_LOST && germlineStatus == HOM_DELETION);
    }
}
