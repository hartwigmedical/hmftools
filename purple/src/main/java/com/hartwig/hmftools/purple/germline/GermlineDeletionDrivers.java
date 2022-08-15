package com.hartwig.hmftools.purple.germline;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.drivercatalog.CNADrivers.MAX_COPY_NUMBER_DEL;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
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
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDetectionMethod;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.GermlineStatus;

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

    public void findDeletions(final List<PurpleCopyNumber> copyNumbers, final List<ObservedRegion> fittedRegions)
    {
        int cnChrStartIndex = 0;
        int cnChrEndIndex = 0;
        String currentCnChromosome = "";

        for(ObservedRegion region : fittedRegions)
        {
            if(region.germlineStatus() != HOM_DELETION && region.germlineStatus() != HET_DELETION)
                continue;

            PurpleCopyNumber matchedCopyNumber = null;

            if(!copyNumbers.isEmpty())
            {
                String chromosome = region.chromosome();
                int regionStart = region.start();
                int regionEnd = region.end();

                // find the overlapping / matching copy number region (skipped for germline-only)
                // both copy numbers and fitted regions are ordered by chromosome and position
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

                    if(positionsOverlap(copyNumber.start(), (int) copyNumber.end(), regionStart, regionEnd))
                    {
                        matchedCopyNumber = copyNumber;
                        break;
                    }
                }
            }

            findOverlappingDriverGene(region, matchedCopyNumber);
        }
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

    private void findOverlappingDriverGene(final ObservedRegion region, final PurpleCopyNumber matchedCopyNumber)
    {
        // now find genes
        List<GeneData> geneDataList = mGeneDataCache.getChrGeneDataMap().get(region.chromosome());

        if(geneDataList == null)
            return;

        int regionLowerPos = region.start() - GERMLINE_DEL_GENE_BUFFER;
        int regionHighPos = region.end() + GERMLINE_DEL_GENE_BUFFER;

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
                DriverGene driverGene = mDriverGenes.stream().filter(y -> y.gene().equals(geneData.GeneName)).findFirst().orElse(null);

                if(driverGene != null)
                {
                    driverGenes.add(driverGene);
                    overlappingGenes.add(geneData);
                }
            }
        }

        if(overlappingGenes.isEmpty() && driverGenes.isEmpty())
            return;

        int cohortFrequency = mCohortFrequency.getRegionFrequency(
                region.chromosome(), region.start(), region.end(), GERMLINE_DEL_REGION_MATCH_BUFFER);

        final List<String> filters = checkFilters(region, matchedCopyNumber, cohortFrequency);

        int exonRankMin = 0;
        int exonRankMax = 0;

        List<GeneData> deletedGenes = Lists.newArrayList();

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

            // take the first if more than one overlapping gene
            if(transcripts.isEmpty())
            {
                exonRankMin = overlappedExons.stream().mapToInt(x -> x.Rank).min().orElse(0);
                exonRankMax = overlappedExons.stream().mapToInt(x -> x.Rank).max().orElse(0);
            }

            transcripts.add(transData);
            deletedGenes.add(geneData);
        }

        if(transcripts.isEmpty())
            return;

        double germlineCopyNumber = region.observedNormalRatio() * 2;
        double tumorCopyNumber = region.refNormalisedCopyNumber();
        GermlineStatus tumorStatus = tumorCopyNumber < 0.5 ? HOM_DELETION : HET_DELETION;

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

        boolean anyReported = filters.isEmpty() && driverGenes.stream().anyMatch(x -> x.reportGermlineDisruption());

        for(GeneData geneData : deletedGenes)
        {
            mDeletions.add(new GermlineDeletion(
                    geneData.GeneName, region.chromosome(), geneData.KaryotypeBand, region.start(), region.end(),
                    region.depthWindowCount(), exonRankMin, exonRankMax,
                    GermlineDetectionMethod.SEGMENT, region.germlineStatus(), tumorStatus, germlineCopyNumber, tumorCopyNumber,
                    filter, cohortFrequency, anyReported));
        }

        for(TranscriptData transData : transcripts)
        {
            GeneData geneData = overlappingGenes.stream().filter(x -> x.GeneId.equals(transData.GeneId)).findFirst().orElse(null);
            DriverGene driverGene = driverGenes.stream().filter(x -> x.gene().equals(geneData.GeneName)).findFirst().orElse(null);

            boolean reported = filters.isEmpty() && driverGene.reportGermlineDisruption();

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
