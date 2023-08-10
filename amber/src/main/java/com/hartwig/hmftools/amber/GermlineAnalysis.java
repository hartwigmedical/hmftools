package com.hartwig.hmftools.amber;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.BaseDepthFactory.fromAmberSite;
import static com.hartwig.hmftools.common.utils.collection.Multimaps.filterEntries;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SamReaderFactory;

public class GermlineAnalysis
{
    private final AmberConfig mConfig;
    private final HetNormalEvidence mHetNormalEvidence;
    private final ListMultimap<Chromosome, BaseDepth> mSnpCheckedLoci;
    private final ListMultimap<Chromosome, BaseDepth> mHomozygousLoci;
    private final ListMultimap<Chromosome, BaseDepth> mHeterozygousLoci;
    private final List<RegionOfHomozygosity> mRegionsOfHomozygosity;
    private final double mConsanguinityProportion;

    @Nullable private final Chromosome mUniparentalDisomy;

    public GermlineAnalysis(
            final AmberConfig config, SamReaderFactory readerFactory, ListMultimap<Chromosome,AmberSite> chromosomeSites)
            throws InterruptedException, IOException
    {
        mConfig = config;

        final Predicate<BaseDepth> isValidFilter = BaseDepth::isValid;
        Predicate<BaseDepth> homozygousFilter = new NormalHomozygousFilter().and(isValidFilter);
        Predicate<BaseDepth> heterozygousFilter = new NormalHeterozygousFilter(mConfig.MinHetAfPercent, mConfig.MaxHetAfPercent).and(isValidFilter);
        Predicate<BaseDepth> snpCheckFilter = new SnpCheckFilter(chromosomeSites);

        mHetNormalEvidence = new HetNormalEvidence();

        // Primary Reference Data
        ListMultimap<Chromosome, BaseDepth> unfilteredLoci = germlineDepth(readerFactory, mConfig.ReferenceBams.get(0), chromosomeSites);

        final Predicate<BaseDepth> depthFilter = new BaseDepthFilter(mConfig.MinDepthPercent, mConfig.MaxDepthPercent, unfilteredLoci);
        mSnpCheckedLoci = filterEntries(unfilteredLoci, snpCheckFilter);
        mHomozygousLoci = filterEntries(unfilteredLoci, depthFilter.and(homozygousFilter));
        var primaryHeterozygousLoci = filterEntries(unfilteredLoci, depthFilter.and(heterozygousFilter));
        mHetNormalEvidence.add(mConfig.primaryReference(), primaryHeterozygousLoci.values());

        // Additional Reference Data
        for(int i = 1; i < mConfig.ReferenceIds.size(); i++)
        {
            final String sample = mConfig.ReferenceIds.get(i);
            final String sampleBam = mConfig.ReferenceBams.get(i);
            final Collection<BaseDepth> additional = germlineDepth(readerFactory, sampleBam, mHetNormalEvidence.intersection()).values();
            final Predicate<BaseDepth> filter = new BaseDepthFilter(mConfig.MinDepthPercent, mConfig.MaxDepthPercent, additional);
            final Collection<BaseDepth> additionalHetNormal = additional.stream().filter(filter.and(heterozygousFilter)).collect(toList());
            mHetNormalEvidence.add(sample, additionalHetNormal);
        }

        if(mConfig.WriteUnfilteredGermline)
            mHeterozygousLoci = unfilteredLoci;
        else
            mHeterozygousLoci = filterEntries(primaryHeterozygousLoci, mHetNormalEvidence.intersectionFilter());

        AMB_LOGGER.info("{} heterozygous, {} homozygous in reference bams", mHeterozygousLoci.size(), mHomozygousLoci.size());

        RegionOfHomozygosityFinder rohFinder = new RegionOfHomozygosityFinder(mConfig.RefGenVersion, mConfig.MinDepthPercent, mConfig.MaxDepthPercent);
        mRegionsOfHomozygosity = rohFinder.findRegions(unfilteredLoci);

        mConsanguinityProportion = ConsanguinityAnalyser.calcConsanguinityProportion(mRegionsOfHomozygosity);
        mUniparentalDisomy = ConsanguinityAnalyser.findUniparentalDisomy(mRegionsOfHomozygosity);
    }

    public ListMultimap<Chromosome, BaseDepth> getSnpCheckedLoci() { return mSnpCheckedLoci; }
    public ListMultimap<Chromosome, BaseDepth> getHomozygousLoci() { return mHomozygousLoci; }
    public ListMultimap<Chromosome, BaseDepth> getHeterozygousLoci() { return mHeterozygousLoci; }
    public List<RegionOfHomozygosity> getRegionsOfHomozygosity() { return mRegionsOfHomozygosity; }
    public double getConsanguinityProportion() { return mConsanguinityProportion; }

    @Nullable
    Chromosome getUniparentalDisomy() { return mUniparentalDisomy; }

    private ListMultimap<Chromosome, BaseDepth> germlineDepth(
            final SamReaderFactory readerFactory, final String bamPath,
            final ListMultimap<Chromosome,AmberSite> bedRegionsSortedSet) throws InterruptedException
    {
        AMB_LOGGER.info("processing {} Amber sites in reference bam({})", bedRegionsSortedSet.values().size(), bamPath);

        Map<Chromosome,List<BaseDepth>> chrBaseDepth = Maps.newHashMap();

        for(Map.Entry<Chromosome,AmberSite> entry : bedRegionsSortedSet.entries())
        {
            Chromosome chromosome = entry.getKey();

            List<BaseDepth> positions = chrBaseDepth.get(chromosome);

            if(positions == null)
            {
                positions = Lists.newArrayList();
                chrBaseDepth.put(chromosome, positions);
            }

            positions.add(fromAmberSite(entry.getValue()));
        }


        BamEvidenceReader bamEvidenceReader = new BamEvidenceReader(mConfig);
        bamEvidenceReader.processBam(bamPath, readerFactory, chrBaseDepth);

        ListMultimap<Chromosome,BaseDepth> normalEvidence = ArrayListMultimap.create();

        for(Map.Entry<Chromosome,List<BaseDepth>> entry : chrBaseDepth.entrySet())
        {
            Chromosome chromosome = entry.getKey();
            List<BaseDepth> positions = entry.getValue();

            positions.forEach(x -> normalEvidence.put(chromosome, x));
        }

        /*
        final List<BaseDepth> baseDepths = bedRegionsSortedSet.values().stream().map(BaseDepthFactory::fromAmberSite).collect(Collectors.toList());

        BaseDepthFactory bafFactory = new BaseDepthFactory(mConfig.MinBaseQuality);

        AsyncBamLociReader.processBam(
                bamPath, readerFactory, baseDepths, bafFactory::addEvidence, mConfig.Threads, mConfig.MinMappingQuality);


        final ListMultimap<Chromosome,BaseDepth> normalEvidence = ArrayListMultimap.create();
        baseDepths.forEach(x -> normalEvidence.put(HumanChromosome.fromString(x.chromosome()), x));
        */

        return normalEvidence;
    }
}
