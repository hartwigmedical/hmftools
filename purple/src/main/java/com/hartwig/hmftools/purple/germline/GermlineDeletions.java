package com.hartwig.hmftools.purple.germline;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.drivercatalog.DeletionDrivers.MAX_COPY_NUMBER_DEL;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.ANY;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.VARIANT_NOT_LOST;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting.WILDTYPE_LOST;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_DEL_CN_CONSISTENCY_MACN_PERC;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_DEL_CN_CONSISTENCY_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_DEL_COHORT_FREQ;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_DEL_GENE_BUFFER;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_DEL_NORMAL_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_DEL_REGION_MATCH_BUFFER;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_DEL_REGION_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.WINDOW_SIZE;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDetectionMethod;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.jetbrains.annotations.Nullable;

public class GermlineDeletions
{
    private final List<DriverGene> mDriverGenes;
    private final EnsemblDataCache mGeneDataCache;
    private final GermlineDeletionFrequency mCohortFrequency;

    private final List<GermlineDeletion> mDeletions;
    private final List<DriverCatalog> mDrivers;

    public GermlineDeletions(
            final List<DriverGene> driverGenes, final EnsemblDataCache geneDataCache, final GermlineDeletionFrequency cohortFrequency)
    {
        mGeneDataCache = geneDataCache;
        mDriverGenes = driverGenes;
        mCohortFrequency = cohortFrequency;

        mDeletions = Lists.newArrayList();
        mDrivers = Lists.newArrayList();
    }

    public List<GermlineDeletion> getDeletions() { return mDeletions; }
    public List<DriverCatalog> getDrivers() { return mDrivers; }

    public void findDeletions(
            final List<PurpleCopyNumber> copyNumbers, final List<ObservedRegion> fittedRegions, final List<StructuralVariant> germlineSVs)
    {
        CopyNumberSearchState copyNumberSearchState = new CopyNumberSearchState();

        for(int i = 0; i < fittedRegions.size(); ++i)
        {
            ObservedRegion region = fittedRegions.get(i);

            if(region.germlineStatus() != HOM_DELETION && region.germlineStatus() != HET_DELETION)
                continue;

            PurpleCopyNumber matchedCopyNumber = findMatchingCopyNumber(region, copyNumbers, copyNumberSearchState);

            ObservedRegion nextRegion = i < fittedRegions.size() - 1 ? fittedRegions.get(i + 1) : null;
            if(nextRegion != null && !nextRegion.chromosome().equals(region.chromosome()))
                nextRegion = null;

            MatchedStructuralVariant[] matchingSVs = findMatchingGermlineSVs(region, nextRegion, germlineSVs);

            findOverlappingDriverGene(region, matchedCopyNumber, matchingSVs);
        }
    }

    private class MatchedStructuralVariant
    {
        public final StructuralVariant Variant;
        public final boolean IsStart;

        public MatchedStructuralVariant(final StructuralVariant variant, final boolean isStart)
        {
            Variant = variant;
            IsStart = isStart;
        }

        public int position() { return Variant.position(IsStart); }
    }

    public MatchedStructuralVariant[] findMatchingGermlineSVs(
            final ObservedRegion region, @Nullable final ObservedRegion nextRegion, final List<StructuralVariant> germlineSVs)
    {
        MatchedStructuralVariant[] matchingSVs = new MatchedStructuralVariant[] { null, null };

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int regionStart;
            int regionEnd;
            byte requiredOrientation = (se == SE_START) ? POS_ORIENT : NEG_ORIENT;

            if(se == SE_START)
            {
                regionStart = region.minStart();
                regionEnd = region.maxStart();
            }
            else
            {
                if(nextRegion != null)
                {
                    regionStart = nextRegion.minStart() - 1;
                    regionEnd = nextRegion.maxStart() - 1;
                }
                else
                {
                    regionStart = region.end();
                    regionEnd = region.end();
                }
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
                        continue;

                    if(!legMatchesRegion(leg, region.chromosome(), regionStart, regionEnd, requiredOrientation))
                        continue;

                    if(matchingSVs[se] == null
                    || (matchingSVs[se] != null && selectNewBySvType(variant.type(), matchingSVs[se].Variant.type())))
                    {
                        matchingSVs[se] = new MatchedStructuralVariant(variant, isStart);
                    }
                }
            }
        }

        return matchingSVs;
    }

    private static boolean legMatchesRegion(
            final StructuralVariantLeg leg, final String regionChr, int regionStart, int regionEnd,  byte requiredOrientation)
    {
        return leg.chromosome().equals(regionChr) && positionWithin(leg.position(), regionStart, regionEnd)
                && leg.orientation() == requiredOrientation;
    }

    private static boolean selectNewBySvType(final StructuralVariantType newType, final StructuralVariantType existingType)
    {
        if(existingType == DEL)
            return false;

        return newType == DEL;
    }

    private class CopyNumberSearchState
    {
        public int ChrStartIndex = 0;
        public int ChrEndIndex = 0;
        public String CurrentChromosome = "";

        public CopyNumberSearchState() {}
    }

    private PurpleCopyNumber findMatchingCopyNumber(
            final ObservedRegion region, final List<PurpleCopyNumber> copyNumbers, final CopyNumberSearchState searchState)
    {
        if(copyNumbers.isEmpty())
            return null;

        PurpleCopyNumber matchedCopyNumber = null;

        String chromosome = region.chromosome();
        int regionStart = region.start();
        int regionEnd = region.end();

        // find the overlapping / matching copy number region (skipped for germline-only)
        // both copy numbers and fitted regions are ordered by chromosome and position
        if(!searchState.CurrentChromosome.equals(chromosome))
        {
            for(int index = searchState.ChrStartIndex; index < copyNumbers.size(); ++index)
            {
                PurpleCopyNumber copyNumber = copyNumbers.get(index);

                if(copyNumber.chromosome().equals(chromosome) && !searchState.CurrentChromosome.equals(chromosome))
                {
                    searchState.CurrentChromosome = chromosome;
                    searchState.ChrStartIndex = index;
                }
                else if(searchState.CurrentChromosome.equals(chromosome) && !copyNumber.chromosome().equals(chromosome))
                {
                    searchState.ChrEndIndex = index - 1;
                    break;
                }
            }
        }

        if(searchState.ChrEndIndex < searchState.ChrStartIndex)
            searchState.ChrEndIndex = copyNumbers.size() - 1;

        for(int index = searchState.ChrStartIndex; index <= searchState.ChrEndIndex; ++index)
        {
            PurpleCopyNumber copyNumber = copyNumbers.get(index);

            if(!copyNumber.chromosome().equals(chromosome))
                break;

            if(positionsOverlap(copyNumber.start(), copyNumber.end(), regionStart, regionEnd))
            {
                matchedCopyNumber = copyNumber;
                break;
            }
        }

        return matchedCopyNumber;
    }

    private static final String FILTER_CN_INCONSISTENCY = "INCONSISTENT_CN";
    private static final String FILTER_COHORT_FREQ = "COHORT_FREQ";
    private static final String FILTER_REGION_LENGTH = "MIN_LENGTH";

    private List<String> checkFilters(final ObservedRegion region, final PurpleCopyNumber copyNumber, int cohortFrequency)
    {
        final List<String> filters = Lists.newArrayList();

        if(region.end() - region.start() <= GERMLINE_DEL_REGION_MIN)
        {
            filters.add(FILTER_REGION_LENGTH);
        }

        if(copyNumber != null)
        {
            double cnInconsistency = region.refNormalisedCopyNumber() - copyNumber.majorAlleleCopyNumber();
            double cnLimit =
                    max(GERMLINE_DEL_CN_CONSISTENCY_MIN, copyNumber.majorAlleleCopyNumber() * GERMLINE_DEL_CN_CONSISTENCY_MACN_PERC);

            if(cnInconsistency > cnLimit)
            {
                filters.add(FILTER_CN_INCONSISTENCY);
            }
        }

        if(!filters.contains(FILTER_CN_INCONSISTENCY) && region.observedNormalRatio() > GERMLINE_DEL_NORMAL_RATIO)
        {
            filters.add(FILTER_CN_INCONSISTENCY);
        }

        if(cohortFrequency >= GERMLINE_DEL_COHORT_FREQ)
        {
            filters.add(FILTER_COHORT_FREQ);
        }

        return filters;
    }

    private void findOverlappingDriverGene(
            final ObservedRegion region, final PurpleCopyNumber matchedCopyNumber, final MatchedStructuralVariant[] matchingSVs)
    {
        // now find genes
        List<GeneData> geneDataList = mGeneDataCache.getChrGeneDataMap().get(region.chromosome());

        if(geneDataList == null)
            return;

        int adjustPosStart = matchingSVs[SE_START] != null ? matchingSVs[SE_START].position() : region.start();
        int adjustPosEnd = matchingSVs[SE_END] != null ? matchingSVs[SE_END].position() : region.end();

        int regionLowerPos = adjustPosStart - GERMLINE_DEL_GENE_BUFFER;
        int regionHighPos = adjustPosEnd + GERMLINE_DEL_GENE_BUFFER;

        double tumorCopyNumber = region.refNormalisedCopyNumber();
        GermlineStatus tumorStatus = tumorCopyNumber < 0.5 ? HOM_DELETION : HET_DELETION;

        // find overlapping driver genes
        List<GeneData> overlappingGenes = Lists.newArrayList();
        List<DriverGene> driverGenes = Lists.newArrayList();
        List<TranscriptData> transcripts = Lists.newArrayList();

        for(GeneData geneData : geneDataList)
        {
            if(regionLowerPos > geneData.GeneEnd)
                continue;

            if(regionHighPos < geneData.GeneStart)
                break;

            if(positionsOverlap(geneData.GeneStart, geneData.GeneEnd, regionLowerPos, regionHighPos))
            {
                overlappingGenes.add(geneData);

                DriverGene driverGene = mDriverGenes.stream().filter(y -> y.gene().equals(geneData.GeneName)).findFirst().orElse(null);

                if(driverGene != null)
                {
                    // check requirements on the germline disruption field: WILDTYPE_LOST - requires an LOH for the deletion to be reportable
                    if(driverGene.reportGermlineDeletion() == ANY || driverGene.reportGermlineDeletion() == VARIANT_NOT_LOST)
                        driverGenes.add(driverGene);
                    else if(driverGene.reportGermlineDeletion() == WILDTYPE_LOST && tumorStatus == HOM_DELETION)
                        driverGenes.add(driverGene);
                }
            }
        }

        if(overlappingGenes.isEmpty())
            return;

        int cohortFrequency = mCohortFrequency.getRegionFrequency(
                region.chromosome(), region.start(), region.end(), GERMLINE_DEL_REGION_MATCH_BUFFER);

        final List<String> filters = checkFilters(region, matchedCopyNumber, cohortFrequency);

        List<GeneData> deletedGenes = Lists.newArrayList();
        List<int[]> deletedExonRanges = Lists.newArrayList();

        for(GeneData geneData : overlappingGenes)
        {
            TranscriptData transData = mGeneDataCache.getTranscriptData(geneData.GeneId, "");

            if(transData == null)
                continue;

            List<ExonData> overlappedExons = transData.exons().stream()
                    .filter(x -> positionsOverlap(x.Start, x.End, regionLowerPos, regionHighPos))
                    .collect(Collectors.toList());

            if(overlappedExons.isEmpty())
                continue;

            PPL_LOGGER.trace("region({}: {}-{}) overlaps gene({}) exons({})",
                    region.chromosome(), region.start(), region.end(), geneData.GeneName, overlappedExons.size());

            int exonRankMin = overlappedExons.stream().mapToInt(x -> x.Rank).min().orElse(0);
            int exonRankMax = overlappedExons.stream().mapToInt(x -> x.Rank).max().orElse(0);

            transcripts.add(transData);
            deletedGenes.add(geneData);
            deletedExonRanges.add(new int[] {exonRankMin, exonRankMax} );
        }

        if(transcripts.isEmpty())
            return;

        double germlineCopyNumber = region.observedNormalRatio() * 2;

        String filter;

        if(filters.isEmpty())
        {
            filter = PASS;
        }
        else
        {
            StringJoiner sj = new StringJoiner(";");
            filters.forEach(x -> sj.add(x));
            filter = sj.toString();
        }

        for(int i = 0; i < deletedGenes.size(); ++i)
        {
            GeneData geneData = deletedGenes.get(i);
            int[] deletedExonRange = deletedExonRanges.get(i);

            boolean reportedGene = filters.isEmpty() && driverGenes.stream().anyMatch(x -> x.gene().equals(geneData.GeneName));

            mDeletions.add(new GermlineDeletion(
                    geneData.GeneName, region.chromosome(), geneData.KaryotypeBand, adjustPosStart, adjustPosEnd,
                    region.depthWindowCount(), deletedExonRange[0], deletedExonRange[1],
                    GermlineDetectionMethod.SEGMENT, region.germlineStatus(), tumorStatus, germlineCopyNumber, tumorCopyNumber,
                    filter, cohortFrequency, reportedGene));
        }

        for(TranscriptData transData : transcripts)
        {
            GeneData geneData = overlappingGenes.stream().filter(x -> x.GeneId.equals(transData.GeneId)).findFirst().orElse(null);
            DriverGene driverGene = driverGenes.stream().filter(x -> x.gene().equals(geneData.GeneName)).findFirst().orElse(null);

            boolean reported = filters.isEmpty() && driverGene != null;

            if(!reported)
                continue;

            // only create one record even if multiple sections of the gene are deleted
            if(mDrivers.stream().anyMatch(x -> x.gene().equals(geneData.GeneName) && x.transcript().equals(transData.TransName)))
                continue;

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
}
