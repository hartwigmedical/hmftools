package com.hartwig.hmftools.amber;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.PositionEvidenceChecker.fromAmberSite;
import static com.hartwig.hmftools.common.utils.Multimaps.filterEntries;

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
    private final ListMultimap<Chromosome, PositionEvidence> mSnpCheckedLoci;
    private final ListMultimap<Chromosome, PositionEvidence> mHomozygousLoci;
    private final ListMultimap<Chromosome, PositionEvidence> mHeterozygousLoci;
    private final List<RegionOfHomozygosity> mRegionsOfHomozygosity;
    private final double mConsanguinityProportion;

    @Nullable private final Chromosome mUniparentalDisomy;

    public GermlineAnalysis(
            final AmberConfig config, SamReaderFactory readerFactory, ListMultimap<Chromosome,AmberSite> chrAmberSites)
            throws InterruptedException, IOException
    {
        mConfig = config;

        final Predicate<PositionEvidence> isValidFilter = PositionEvidence::isValid;
        Predicate<PositionEvidence> homozygousFilter = new NormalHomozygousFilter().and(isValidFilter);
        Predicate<PositionEvidence> heterozygousFilter = new NormalHeterozygousFilter(mConfig.MinHetAfPercent, mConfig.MaxHetAfPercent).and(isValidFilter);
        Predicate<PositionEvidence> snpCheckFilter = new SnpCheckFilter(chrAmberSites);

        mHetNormalEvidence = new HetNormalEvidence();

        // Primary Reference Data
        ListMultimap<Chromosome, PositionEvidence> unfilteredLoci = germlineDepth(readerFactory, mConfig.ReferenceBams.get(0), chrAmberSites);

        final Predicate<PositionEvidence> depthFilter = new BaseDepthFilter(mConfig.MinDepthPercent, mConfig.MaxDepthPercent, unfilteredLoci);
        mSnpCheckedLoci = filterEntries(unfilteredLoci, snpCheckFilter);
        mHomozygousLoci = filterEntries(unfilteredLoci, depthFilter.and(homozygousFilter));
        var primaryHeterozygousLoci = filterEntries(unfilteredLoci, depthFilter.and(heterozygousFilter));
        mHetNormalEvidence.add(mConfig.primaryReference(), primaryHeterozygousLoci.values());

        // Additional Reference Data
        for(int i = 1; i < mConfig.ReferenceIds.size(); i++)
        {
            final String sample = mConfig.ReferenceIds.get(i);
            final String sampleBam = mConfig.ReferenceBams.get(i);
            final Collection<PositionEvidence> additional = germlineDepth(readerFactory, sampleBam, mHetNormalEvidence.intersection()).values();
            final Predicate<PositionEvidence> filter = new BaseDepthFilter(mConfig.MinDepthPercent, mConfig.MaxDepthPercent, additional);
            final Collection<PositionEvidence> additionalHetNormal = additional.stream().filter(filter.and(heterozygousFilter)).collect(toList());
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

    public ListMultimap<Chromosome,PositionEvidence> getSnpCheckedLoci() { return mSnpCheckedLoci; }
    public ListMultimap<Chromosome,PositionEvidence> getHomozygousLoci() { return mHomozygousLoci; }
    public ListMultimap<Chromosome,PositionEvidence> getHeterozygousLoci() { return mHeterozygousLoci; }

    public List<RegionOfHomozygosity> getRegionsOfHomozygosity() { return mRegionsOfHomozygosity; }
    public double getConsanguinityProportion() { return mConsanguinityProportion; }

    @Nullable
    Chromosome getUniparentalDisomy() { return mUniparentalDisomy; }

    private ListMultimap<Chromosome, PositionEvidence> germlineDepth(
            final SamReaderFactory readerFactory, final String bamPath,
            final ListMultimap<Chromosome,AmberSite> chrAmberSites) throws InterruptedException
    {
        AMB_LOGGER.info("processing {} Amber sites in reference bam({})", chrAmberSites.values().size(), bamPath);

        Map<Chromosome,List<PositionEvidence>> chrPositionEvidence = Maps.newHashMap();

        for(Map.Entry<Chromosome,AmberSite> entry : chrAmberSites.entries())
        {
            Chromosome chromosome = entry.getKey();

            List<PositionEvidence> positions = chrPositionEvidence.get(chromosome);

            if(positions == null)
            {
                positions = Lists.newArrayList();
                chrPositionEvidence.put(chromosome, positions);
            }

            positions.add(fromAmberSite(entry.getValue()));
        }

        BamEvidenceReader bamEvidenceReader = new BamEvidenceReader(mConfig);
        bamEvidenceReader.processBam(bamPath, readerFactory, chrPositionEvidence);

        ListMultimap<Chromosome, PositionEvidence> normalEvidence = ArrayListMultimap.create();

        for(Map.Entry<Chromosome,List<PositionEvidence>> entry : chrPositionEvidence.entrySet())
        {
            Chromosome chromosome = entry.getKey();
            List<PositionEvidence> positions = entry.getValue();

            positions.forEach(x -> normalEvidence.put(chromosome, x));
        }

        return normalEvidence;
    }
}
