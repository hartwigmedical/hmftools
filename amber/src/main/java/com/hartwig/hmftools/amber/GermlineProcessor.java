package com.hartwig.hmftools.amber;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.common.utils.collection.Multimaps.filterEntries;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.function.Predicate;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.NormalHeterozygousFilter;
import com.hartwig.hmftools.common.amber.NormalHomozygousFilter;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SamReaderFactory;

public class GermlineProcessor
{
    private final AmberConfig mConfig;
    private final ExecutorService mExecutorService;
    private final AmberHetNormalEvidence mHetNormalEvidence;
    private final ListMultimap<Chromosome, BaseDepth> mSnpCheckedLoci;
    private final ListMultimap<Chromosome, BaseDepth> mHomozygousLoci;
    private final ListMultimap<Chromosome, BaseDepth> mHeterozygousLoci;
    private final List<RegionOfHomozygosity> mRegionsOfHomozygosity;
    private final double mConsanguinityProportion;
    @Nullable private final Chromosome mUniparentalDisomy;

    AmberHetNormalEvidence getHetEvidence() { return mHetNormalEvidence; }

    ListMultimap<Chromosome, BaseDepth> getSnpCheckedLoci() { return mSnpCheckedLoci; }
    ListMultimap<Chromosome, BaseDepth> getHomozygousLoci() { return mHomozygousLoci; }
    ListMultimap<Chromosome, BaseDepth> getHeterozygousLoci() { return mHeterozygousLoci; }

    List<RegionOfHomozygosity> getRegionsOfHomozygosity() { return mRegionsOfHomozygosity; }
    double getConsanguinityProportion() { return mConsanguinityProportion; }

    @Nullable
    Chromosome getUniparentalDisomy() { return mUniparentalDisomy; }

    GermlineProcessor(final AmberConfig config, SamReaderFactory readerFactory, ExecutorService executorService, ListMultimap<Chromosome,AmberSite> chromosomeSites)
            throws InterruptedException, ExecutionException
    {
        mConfig = config;
        mExecutorService = executorService;

        final Predicate<BaseDepth> isValidFilter = BaseDepth::isValid;
        Predicate<BaseDepth> homozygousFilter = new NormalHomozygousFilter().and(isValidFilter);
        Predicate<BaseDepth> heterozygousFilter = new NormalHeterozygousFilter(mConfig.MinHetAfPercent, mConfig.MaxHetAfPercent).and(isValidFilter);
        Predicate<BaseDepth> snpCheckFilter = new SnpCheckFilter(chromosomeSites);

        mHetNormalEvidence = new AmberHetNormalEvidence();

        // Primary Reference Data
        final ListMultimap<Chromosome, BaseDepth> unfilteredLoci = normalDepth(readerFactory, mConfig.ReferenceBamPath.get(0),
                chromosomeSites);
        final Predicate<BaseDepth> depthFilter = new BaseDepthFilter(mConfig.MinDepthPercent, mConfig.MaxDepthPercent, unfilteredLoci);
        mSnpCheckedLoci = filterEntries(unfilteredLoci, snpCheckFilter);
        mHomozygousLoci = filterEntries(unfilteredLoci, depthFilter.and(homozygousFilter));
        var primaryHeterozygousLoci = filterEntries(unfilteredLoci, depthFilter.and(heterozygousFilter));
        mHetNormalEvidence.add(mConfig.primaryReference(), primaryHeterozygousLoci.values());

        // Additional Reference Data
        for(int i = 1; i < mConfig.ReferenceIds.size(); i++)
        {
            final String sample = mConfig.ReferenceIds.get(i);
            final String sampleBam = mConfig.ReferenceBamPath.get(i);
            final Collection<BaseDepth> additional = normalDepth(readerFactory, sampleBam, mHetNormalEvidence.intersection()).values();
            final Predicate<BaseDepth> filter = new BaseDepthFilter(mConfig.MinDepthPercent, mConfig.MaxDepthPercent, additional);
            final Collection<BaseDepth> additionalHetNormal = additional.stream().filter(filter.and(heterozygousFilter)).collect(toList());
            mHetNormalEvidence.add(sample, additionalHetNormal);
        }

        mHeterozygousLoci = filterEntries(primaryHeterozygousLoci, mHetNormalEvidence.intersectionFilter());

        AMB_LOGGER.info("{} loci in reference bams", mHeterozygousLoci.size());

        RegionOfHomozygosityFinder rohFinder = new RegionOfHomozygosityFinder(mConfig.refGenomeVersion, mConfig.MinDepthPercent, mConfig.MaxDepthPercent);
        mRegionsOfHomozygosity = rohFinder.findRegions(unfilteredLoci);

        mConsanguinityProportion = ConsanguinityAnalyser.calcConsanguinityProportion(mRegionsOfHomozygosity);
        mUniparentalDisomy = ConsanguinityAnalyser.findUniparentalDisomy(mRegionsOfHomozygosity);
    }

    private ListMultimap<Chromosome, BaseDepth> normalDepth(
            final SamReaderFactory readerFactory, final String bamPath,
            final ListMultimap<Chromosome, AmberSite> bedRegionsSortedSet) throws InterruptedException, ExecutionException
    {
        final int partitionSize = Math.max(mConfig.minPartition(), bedRegionsSortedSet.size() / mConfig.ThreadCount);

        AMB_LOGGER.info("Processing {} potential sites in reference bam {}", bedRegionsSortedSet.values().size(), bamPath);
        final AmberTaskCompletion completion = new AmberTaskCompletion();

        final List<Future<BaseDepthEvidence>> futures = new ArrayList<>();
        for(final Chromosome contig : bedRegionsSortedSet.keySet())
        {
            for(final List<AmberSite> inner : Lists.partition(new ArrayList<>(bedRegionsSortedSet.get(contig)), partitionSize))
            {
                final BaseDepthEvidence evidence = new BaseDepthEvidence(mConfig.typicalReadDepth(),
                        mConfig.MinMappingQuality, mConfig.MinBaseQuality,
                        inner.get(0).chromosome(), bamPath, readerFactory, inner);

                futures.add(mExecutorService.submit(completion.task(evidence)));
            }
        }

        final ListMultimap<Chromosome, BaseDepth> normalEvidence = ArrayListMultimap.create();
        AmberUtils.getFuture(futures).forEach(x -> normalEvidence.putAll(HumanChromosome.fromString(x.contig()), x.evidence()));
        return normalEvidence;
    }
}
