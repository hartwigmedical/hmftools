package com.hartwig.hmftools.purple.drivers;

import static java.lang.Math.copySign;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.purple.region.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.UNKNOWN;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.GERMLINE_DEL_CN_CONSISTENCY_MACN_PERC;
import static com.hartwig.hmftools.purple.config.PurpleConstants.GERMLINE_DEL_CN_CONSISTENCY_MIN;
import static com.hartwig.hmftools.purple.config.PurpleConstants.GERMLINE_DEL_COHORT_FREQ;
import static com.hartwig.hmftools.purple.config.PurpleConstants.GERMLINE_DEL_GENE_BUFFER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.GERMLINE_DEL_NORMAL_RATIO;
import static com.hartwig.hmftools.purple.config.PurpleConstants.GERMLINE_DEL_REGION_MATCH_BUFFER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.GERMLINE_DEL_REGION_MIN;

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
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GermlineDeletion;
import com.hartwig.hmftools.common.purple.gene.GermlineDetectionMethod;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;

public class GermlineDeletionDrivers
{
    private final List<DriverGene> mDriverGenes;
    private final EnsemblDataCache mGeneDataCache;
    private final GermlineDeletionFrequency mCohortFrequency;

    private final List<GermlineDeletion> mDeletions;
    private final List<DriverCatalog> mDrivers;

    public GermlineDeletionDrivers(
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

    public void findDeletions(final List<PurpleCopyNumber> copyNumbers, final List<FittedRegion> fittedRegions)
    {
        int cnChrStartIndex = 0;
        int cnChrEndIndex = 0;
        String currentCnChromosome = "";

        for(FittedRegion region : fittedRegions)
        {
            if(region.germlineStatus() != HOM_DELETION && region.germlineStatus() != HET_DELETION)
                continue;

            String chromosome = region.chromosome();
            int regionStart = (int)region.start();
            int regionEnd = (int)region.end();

            // find the overlapping / matching copy number region
            // both copy numbers and fitted regions are ordered by chromosome and position
            PurpleCopyNumber matchedCopyNumber = null;

            if(!currentCnChromosome.equals(chromosome))
            {
                for(int index = cnChrStartIndex; index < copyNumbers.size(); ++index)
                {
                    PurpleCopyNumber copyNumber = copyNumbers.get(index);

                    if(copyNumber.chromosome().equals(chromosome) && !currentCnChromosome.equals(chromosome))
                    {
                        currentCnChromosome = chromosome;
                        cnChrStartIndex = index;
                    }
                    else if(currentCnChromosome.equals(chromosome) && !copyNumber.chromosome().equals(chromosome))
                    {
                        cnChrEndIndex = index - 1;
                        break;
                    }
                }
            }

            if(cnChrEndIndex < cnChrStartIndex)
                cnChrEndIndex = copyNumbers.size() - 1;

            for(int index = cnChrStartIndex; index <= cnChrEndIndex; ++index)
            {
                PurpleCopyNumber copyNumber = copyNumbers.get(index);

                if(!copyNumber.chromosome().equals(chromosome))
                    break;

                if(positionsOverlap((int)copyNumber.start(), (int)copyNumber.end(), regionStart, regionEnd))
                {
                    matchedCopyNumber = copyNumber;
                    break;
                }
            }

            if(matchedCopyNumber == null)
                continue;

            int cohortFrequency = mCohortFrequency.getRegionFrequency(
                    region.chromosome(), (int)region.start(), (int)region.end(), GERMLINE_DEL_REGION_MATCH_BUFFER);

            final List<String> filters = checkFilters(region, matchedCopyNumber, cohortFrequency);

            findOverlappingDriverGene(region, filters, cohortFrequency);
        }
    }

    private static final String FILTER_CN_INCONSISTENCY = "INCONSISTENT_CN";
    private static final String FILTER_COHORT_FREQ = "COHORT_FREQ";
    private static final String FILTER_REGION_LENGTH = "MIN_LENGTH";

    private List<String> checkFilters(final FittedRegion region, final PurpleCopyNumber copyNumber, int cohortFrequency)
    {
        final List<String> filters = Lists.newArrayList();

        if(region.end() - region.start() <= GERMLINE_DEL_REGION_MIN)
        {
            filters.add(FILTER_REGION_LENGTH);
        }

        // Filter (pass, minLength, inconsistentTumorCN,highNormalRatio)
        double cnInconsistency = region.refNormalisedCopyNumber() - copyNumber.majorAlleleCopyNumber();
        double cnLimit = max(GERMLINE_DEL_CN_CONSISTENCY_MIN, copyNumber.majorAlleleCopyNumber() * GERMLINE_DEL_CN_CONSISTENCY_MACN_PERC);

        if(cnInconsistency > cnLimit)
        {
            filters.add(FILTER_CN_INCONSISTENCY);
        }

        if(region.observedNormalRatio() > GERMLINE_DEL_NORMAL_RATIO)
        {
            filters.add(FILTER_CN_INCONSISTENCY);
        }

        if(cohortFrequency >= GERMLINE_DEL_COHORT_FREQ)
        {
            filters.add(FILTER_COHORT_FREQ);
        }

        return filters;
    }

    private void findOverlappingDriverGene(final FittedRegion region, final List<String> filters, int cohortFrequency)
    {
        // now find genes
        List<GeneData> geneDataList = mGeneDataCache.getChrGeneDataMap().get(region.chromosome());

        if(geneDataList == null)
            return;

        int regionLowerPos = (int)region.start() - GERMLINE_DEL_GENE_BUFFER;
        int regionHighPos = (int)region.end() + GERMLINE_DEL_GENE_BUFFER;

        // find overlapping driver genes
        GeneData overlappingGeneData = null;
        DriverGene driverGene = null;

        for(GeneData geneData : geneDataList)
        {
            if(regionLowerPos > geneData.GeneEnd)
                continue;

            if(regionHighPos < geneData.GeneStart)
                break;

            if(positionsOverlap(geneData.GeneStart, geneData.GeneEnd, regionLowerPos, regionHighPos))
            {
                overlappingGeneData = geneData;
                driverGene = mDriverGenes.stream().filter(y -> y.gene().equals(geneData.GeneName)).findFirst().orElse(null);
                break;
            }
        }

        if(overlappingGeneData == null || driverGene == null)
            return;

        TranscriptData transData = mGeneDataCache.getTranscriptData(overlappingGeneData.GeneId, "");

        if(transData == null)
            return;

        List<ExonData> overlappedExons = transData.exons().stream()
                .filter(x -> positionsOverlap(x.Start, x.End, regionLowerPos, regionHighPos))
                .collect(Collectors.toList());

        if(overlappedExons.isEmpty())
            return;

        PPL_LOGGER.debug("region({}: {}-{}) overlaps gene({}) exons({})",
                region.chromosome(), region.start(), region.end(), overlappingGeneData.GeneName, overlappedExons.size());

        int exonRankMin = overlappedExons.stream().mapToInt(x -> x.Rank).min().orElse(0);
        int exonRankMax = overlappedExons.stream().mapToInt(x -> x.Rank).max().orElse(0);
        double germlineCopyNumber = region.observedNormalRatio() * 2;
        double tumorCopyNumber = region.refNormalisedCopyNumber();
        GermlineStatus tumorStatus = tumorCopyNumber < 0.5 ? HOM_DELETION : HET_DELETION;
        boolean reported = filters.isEmpty() && driverGene.reportGermlineDisruption();

        String filter;

        if(filters.isEmpty())
        {
            filter = "PASS";
        }
        else
        {
            StringJoiner sj = new StringJoiner(";");
            filters.forEach(x -> sj.add(x));
            filter = sj.toString();
        }

        mDeletions.add(new GermlineDeletion(
                overlappingGeneData.GeneName, region.chromosome(), (int)region.start(), (int)region.end(),
                region.depthWindowCount(), exonRankMin, exonRankMax,
                GermlineDetectionMethod.SEGMENT, region.germlineStatus(), tumorStatus, germlineCopyNumber, tumorCopyNumber,
                filter, cohortFrequency, reported));

        // create a driver record
        DriverCatalog driverCatalog = ImmutableDriverCatalog.builder()
                .chromosome(region.chromosome())
                .chromosomeBand(overlappingGeneData.KaryotypeBand)
                .gene(overlappingGeneData.GeneName)
                .transcript(transData.TransName)
                .isCanonical(transData.IsCanonical)
                .driver(DriverType.GERMLINE)
                .category(driverGene.likelihoodType())
                .driverLikelihood(1)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .biallelic(false)
                .minCopyNumber(region.refNormalisedCopyNumber())
                .maxCopyNumber(region.refNormalisedCopyNumber())
                .likelihoodMethod(LikelihoodMethod.GERMLINE)
                .build();

        mDrivers.add(driverCatalog);
    }
}
