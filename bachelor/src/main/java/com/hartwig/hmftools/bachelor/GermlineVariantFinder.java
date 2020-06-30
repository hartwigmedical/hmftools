package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.ExternalDBFilters.isBenign;
import static com.hartwig.hmftools.bachelor.ExternalDBFilters.stripTranscriptVersion;
import static com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant.MATCH_TYPE_NONE;
import static com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant.MATCH_TYPE_REQUIRED_EFFECT;
import static com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant.MATCH_TYPE_WHITELIST;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.bachelor.datamodel.GeneIdentifier;
import com.hartwig.hmftools.bachelor.datamodel.Program;
import com.hartwig.hmftools.bachelor.datamodel.ProgramPanel;
import com.hartwig.hmftools.bachelor.datamodel.SnpEffect;
import com.hartwig.hmftools.bachelor.datamodel.VariantException;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelor.types.VariantFilter;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class GermlineVariantFinder
{
    private String mName;

    private List<String> mRequiredEffects;
    private List<String> mPanelTranscripts;

    private Map<String,List<VariantFilter>> mWhitelistFilters;
    private Map<String,List<VariantFilter>> mBlacklistFilters;

    private Map<String, HmfTranscriptRegion> mAllGenesMap;
    private Map<String, HmfTranscriptRegion> mAllTranscriptsMap;

    private Set<HmfTranscriptRegion> mTranscriptRegions;

    private List<BachelorGermlineVariant> mVariants;

    private static final Logger LOGGER = LogManager.getLogger(GermlineVariantFinder.class);

    GermlineVariantFinder()
    {
        mName = "";
        mRequiredEffects = null;
        mPanelTranscripts = null;

        mWhitelistFilters = Maps.newHashMap();
        mBlacklistFilters = Maps.newHashMap();

        mVariants = Lists.newArrayList();

        mTranscriptRegions = Sets.newHashSet();

        initialiseGeneData();
    }

    List<BachelorGermlineVariant> getVariants() { return mVariants; }

    boolean loadConfig(Map<String, Program> input)
    {
        if(input.values().isEmpty())
            return false;

        final Program program = input.values().iterator().next();

        mName = program.getName();

        if(program.getPanel().isEmpty())
            return false;

        final Multimap<String, String> geneToEnsemblMap = HashMultimap.create();

        program.getPanel()
                .stream()
                .map(ProgramPanel::getGene)
                .flatMap(Collection::stream)
                .forEach(g -> geneToEnsemblMap.put(g.getName(), g.getEnsembl()));

        final ProgramPanel panel = program.getPanel().get(0);

        final List<GeneIdentifier> genes = panel.getGene();

        // take up a collection of the effects to search for
        mRequiredEffects = panel.getSnpEffect().stream().map(SnpEffect::value).collect(Collectors.toList());
        mPanelTranscripts = genes.stream().map(GeneIdentifier::getEnsembl).collect(Collectors.toList());

        // update query targets
        for (final GeneIdentifier g : genes)
        {
            final HmfTranscriptRegion region = mAllTranscriptsMap.get(g.getEnsembl());

            if (region == null)
            {
                final HmfTranscriptRegion namedRegion = mAllGenesMap.get(g.getName());

                if (namedRegion == null)
                {
                    LOGGER.warn("Program {} gene {} non-canonical transcript {} couldn't find region, transcript will be skipped",
                            program.getName(), g.getName(), g.getEnsembl());

                    // just skip this gene for now
                }
                else
                {
                    mTranscriptRegions.add(namedRegion);
                }
            }
            else
            {
                mTranscriptRegions.add(region);
            }
        }

        // merge XML white and black lists into the same format
        if(program.getBlacklist() != null)
        {
            for (VariantException exception : program.getBlacklist().getExclusion())
            {
                VariantFilter filter = loadVariantFilter(exception);

                List<VariantFilter> filters = mBlacklistFilters.get(filter.Gene);

                if(filters == null)
                {
                    filters = Lists.newArrayList();
                    mBlacklistFilters.put(filter.Gene, filters);
                }

                filters.add(filter);
            }
        }

        if(program.getWhitelist() != null)
        {
            for(VariantException exception : program.getWhitelist().getInclusion())
            {
                VariantFilter filter = loadVariantFilter(exception);

                List<VariantFilter> filters = mWhitelistFilters.get(filter.Gene);

                if(filters == null)
                {
                    filters = Lists.newArrayList();
                    mWhitelistFilters.put(filter.Gene, filters);
                }

                filters.add(filter);
            }
        }

        return true;
    }

    private VariantFilter loadVariantFilter(final VariantException variantException)
    {
        String hgvsProtein = variantException.getHGVSP() != null ? variantException.getHGVSP() : "";
        int minCodon = variantException.getMinCodon() != null ? variantException.getMinCodon().intValue() : -1;

        String chromosome = "";
        long position = 0;
        String ref = "";
        String alt = "";

        if(variantException.getPosition() != null)
        {
            String[] varDetails = variantException.getPosition().split(":");

            if(varDetails.length == 4)
            {
                chromosome = varDetails[0];
                position = Long.parseLong(varDetails[1]);
                ref = varDetails[2];
                alt = varDetails[3];
            }
        }

        final String gene = variantException.getGene().getName();

        return new VariantFilter(gene, "", chromosome, position,
                ref, alt, NONSENSE_OR_FRAMESHIFT, hgvsProtein, "", "", minCodon);
    }

    void addExternalFilters(final List<VariantFilter> allFilters)
    {
        // split into white and black list based on the coding effect
        for(VariantFilter filter : allFilters)
        {
            List<VariantFilter> filters;

            if(filter.Effect == NONSENSE_OR_FRAMESHIFT || filter.Effect == CodingEffect.SPLICE)
            {
                filters = mBlacklistFilters.get(filter.Gene);

                if(filters == null)
                {
                    filters = Lists.newArrayList();
                    mBlacklistFilters.put(filter.Gene, filters);
                }
            }
            else
            {
                filters = mWhitelistFilters.get(filter.Gene);

                if(filters == null)
                {
                    filters = Lists.newArrayList();
                    mWhitelistFilters.put(filter.Gene, filters);
                }
            }

            filters.add(filter);

        }
    }

    public String name() { return mName; }

    void processVcfFile(final String sampleId, final VCFFileReader reader, boolean usesIndex)
    {
        mVariants.clear();

        if(usesIndex)
        {
            for (final HmfTranscriptRegion region : mTranscriptRegions)
            {
                final CloseableIterator<VariantContext> query =
                        reader.query(region.chromosome(), (int) region.geneStart(), (int) region.geneEnd());

                while (query.hasNext())
                {
                    final VariantContext variant = query.next();
                    processVariant(variant, sampleId, region);
                }

                query.close();
            }
        }
        else
        {
            for (final VariantContext variant : reader)
            {
                processVariant(variant, sampleId, null);
            }
        }
    }

    private void processVariant(final VariantContext variant, final String sampleId, HmfTranscriptRegion region)
    {
        if (variant.isFiltered())
            return;

        // we will skip when an ALT is not present in the sample
        final Genotype refGenotype = variant.getGenotype(0);

        if (refGenotype == null || !(refGenotype.isHomVar() || refGenotype.isHet()))
        {
            return;
        }

        final List<SnpEffAnnotation> sampleAnnotations = SnpEffAnnotationFactory.fromContext(variant);

        // search the list of annotations for the correct allele and transcript ID to write to the result file

        // check the sub-conditions now - hotspot locations and gene-transcript IDs
        String matchType = MATCH_TYPE_NONE;

        // first check the transcript

        for (int i = 0; i < sampleAnnotations.size(); ++i)
        {
            final SnpEffAnnotation snpEff = sampleAnnotations.get(i);

            if(!snpEff.allele().equals(refGenotype.getAllele(1).getBaseString()))
                continue; // skip alleles for the second sample (ie the tumor)

            if (!snpEff.isTranscriptFeature())
                continue;

            final String transcriptId = stripTranscriptVersion(snpEff.transcript());

            if (region != null)
            {
                if(!region.transcriptID().equals(transcriptId))
                    continue;
            }
            else
            {
                if(!mPanelTranscripts.contains(transcriptId))
                    continue;
            }

            final String gene = snpEff.gene();
            CodingEffect codingEffect = CodingEffect.effect(gene, snpEff.consequences());

            final String varId = variant.getID();
            final String chromosome = variant.getContig();
            final long position = variant.getStart();
            final String ref = variant.getReference().toString().replaceAll("\\*", "");
            final String alt = snpEff.allele().replaceAll("\\*", "");

            final String effects = snpEff.effects();
            final String hgvsProtein = snpEff.hgvsProtein();
            final String hgvsCoding = snpEff.hgvsCoding();

            for (String requiredEffect : mRequiredEffects)
            {
                if (effects.contains(requiredEffect))
                {
                    LOGGER.debug("Match found: gene({} {}) var({}:{}) ref({}) alt({}) on effect({})",
                            gene, transcriptId, chromosome, position, ref, alt, effects);

                    matchType = MATCH_TYPE_REQUIRED_EFFECT;
                    break;
                }
            }

            VariantFilter matchedFilter = null;

            if (matchType.equals(MATCH_TYPE_REQUIRED_EFFECT) && !mBlacklistFilters.isEmpty())
            {
                // for variants matching the required effects, check whether they should be blacklisted
                // for Clinvar entries this is if the variant is Benign, for other filters any match will cause a blacklist
                final List<Integer> proteinPositions = proteinPosition(snpEff);
                int proteinPosition = !proteinPositions.isEmpty() ? proteinPositions.get(0) : -1;

                List<VariantFilter> filters = mBlacklistFilters.get(gene);

                if(filters != null)
                {
                    for(final VariantFilter filter : filters)
                    {
                        if(filter.blacklistMatch(gene, chromosome, position, ref, alt, proteinPosition))
                        {
                            matchedFilter = filter;

                            if(isBenign(filter.ClinvarSignificance) || filter.ClinvarSignificance.isEmpty())
                            {
                                LOGGER.debug("Gene({}) var({}:{}:{}) ref({}) alt({}) protein({}) blacklisted",
                                        gene, varId, chromosome, position, ref, alt, hgvsProtein);

                                return;
                            }
                        }
                    }
                }
            }

            if (matchType.equals(MATCH_TYPE_NONE) && !mWhitelistFilters.isEmpty())
            {
                List<VariantFilter> filters = mWhitelistFilters.get(gene);

                if(filters != null)
                {
                    for (final VariantFilter filter : filters)
                    {
                        if(filter.whitelistMatch(gene, chromosome, position, ref, alt, codingEffect, hgvsProtein))
                        {
                            LOGGER.debug("Match found: gene({} {}) var({}:{}:{}) ref({}) alt({}) hgvsProtein({}) whitelisted",
                                    gene, transcriptId, varId, chromosome, position, ref, alt, hgvsProtein);
                            matchType = MATCH_TYPE_WHITELIST;
                            matchedFilter = filter;
                            break;
                        }
                    }
                }
            }

            if(matchType.equals(MATCH_TYPE_NONE))
                return;

            String annotationsStr = SnpEffAnnotationFactory.rawAnnotations(variant).get(i);

            boolean isHomozygous = refGenotype.isHom();
            int phredScore = refGenotype.getPL().length >= 1 ? refGenotype.getPL()[0] : 0;

            int alleleIndex = 1;
            for(; alleleIndex < variant.getAlleles().size(); ++alleleIndex)
            {
                if(variant.getAlleles().get(alleleIndex).getBaseString().equals(refGenotype.getAllele(1).getBaseString()))
                    break;
            }

            int germlineAltCount = alleleIndex < refGenotype.getAD().length ? refGenotype.getAD()[alleleIndex] : 0;
            int germlineReadDepth = refGenotype.getDP();

            int tumorAltCount = 0;
            int tumorReadDepth = 0;
            boolean hasDepthInfo = false;

            if (variant.getGenotypes().size() >= 2)
            {
                hasDepthInfo = true;
                final Genotype tumorGenotype = variant.getGenotype(1);
                int[] alleleData = tumorGenotype.getAD();
                tumorAltCount = alleleData[1];
                tumorReadDepth = tumorGenotype.getDP();
            }

            final String codonInfo = snpEff.aaPosAndLength();

            BachelorGermlineVariant germlineVariant = new BachelorGermlineVariant(
                    sampleId, mName, variant.getID(),
                    gene, transcriptId, chromosome, position, ref, alt,
                    codingEffect, effects, annotationsStr, hgvsProtein, isHomozygous, phredScore, hgvsCoding,
                    matchType, codonInfo);

            if(matchedFilter != null)
            {
                germlineVariant.setClinvarData(matchedFilter.ClinvarSignificance, matchedFilter.ClinvarSigInfo);
            }

            germlineVariant.setGermlineData(germlineAltCount, germlineReadDepth);

            if(hasDepthInfo)
                germlineVariant.setReadData(tumorAltCount, tumorReadDepth);

            LOGGER.debug("Adding germline variant({} {}:{}) to list of germline variants",
                    germlineVariant.Gene, germlineVariant.Chromosome, germlineVariant.Position);
            mVariants.add(germlineVariant);
        }
    }

    private static List<Integer> proteinPosition(final SnpEffAnnotation annotation)
    {
        return Arrays.stream(annotation.aaPosAndLength().split("/"))
                .filter(s -> !s.isEmpty())
                .map(Integer::parseInt)
                .collect(Collectors.toList());
    }

    private void initialiseGeneData()
    {
        SortedSetMultimap<String, HmfTranscriptRegion> mGenesByChromosomeMap = HmfGenePanelSupplier.allGenesPerChromosomeMap37();

        mAllGenesMap = Maps.newHashMap();
        for (final HmfTranscriptRegion region : mGenesByChromosomeMap.values())
        {
            mAllGenesMap.put(region.gene(), region);
        }

        mAllTranscriptsMap = Maps.newHashMap();
        for (final HmfTranscriptRegion region : mGenesByChromosomeMap.values())
        {
            mAllTranscriptsMap.put(region.transcriptID(), region);
        }
    }

    /* HOTSPOT logic

        // then check the hotspot location
        if(matchType == NONE && !mHotspots.isEmpty())
        {
            for (final HotspotLocation hotspot : mHotspots)
            {
                if (position != hotspot.getPosition().intValue() || !chromosome.equals(hotspot.getChromosome()))
                    continue;

                if (!ref.equals(hotspot.getRef()) || variant.getAlleles().size() < 2 || !alt.equals(hotspot.getAlt()))
                {
                    continue;
                }

                matchType = HOTSPOT_LOCATION;

                LOGGER.debug("match found: gene({} {}) var({}:{}) ref({}) alt({}) on hotspot location)",
                        gene, transcriptId, chromosome, position, ref, alt, varId);
            }
        }

    */
}
