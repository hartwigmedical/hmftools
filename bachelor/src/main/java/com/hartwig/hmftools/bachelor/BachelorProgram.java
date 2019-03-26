package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.NONE;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.REQUIRED_EFFECT;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.WHITELIST;
import static com.hartwig.hmftools.bachelor.ExternalDBFilters.isBenign;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;

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
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import nl.hartwigmedicalfoundation.bachelor.GeneIdentifier;
import nl.hartwigmedicalfoundation.bachelor.HotspotLocation;
import nl.hartwigmedicalfoundation.bachelor.Program;
import nl.hartwigmedicalfoundation.bachelor.ProgramBlacklist;
import nl.hartwigmedicalfoundation.bachelor.ProgramWhitelist;
import nl.hartwigmedicalfoundation.bachelor.ProgramPanel;
import nl.hartwigmedicalfoundation.bachelor.SnpEffect;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class BachelorProgram
{
    private String mName;

    private List<String> mRequiredEffects;
    private List<String> mPanelTranscripts;

    private Map<String,List<VariantFilter>> mWhitelistFilters;
    private Map<String,List<VariantFilter>> mBlacklistFilters;

    private SortedSetMultimap<String, HmfTranscriptRegion> mGenesByChromosomeMap;
    private Map<String, HmfTranscriptRegion> mAllGenesMap;
    private Map<String, HmfTranscriptRegion> mAllTranscriptsMap;

    private Set<HmfTranscriptRegion> mTranscriptRegions;

    private List<EligibilityReport> mReports;

    private static final Logger LOGGER = LogManager.getLogger(BachelorProgram.class);

    public BachelorProgram()
    {
        mName = "";
        mRequiredEffects = null;
        mPanelTranscripts = null;

        mWhitelistFilters = Maps.newHashMap();
        mBlacklistFilters = Maps.newHashMap();

        mReports = Lists.newArrayList();

        mTranscriptRegions = Sets.newHashSet();

        initialiseGeneData();
    }

    public boolean loadConfig(Map<String, Program> input)
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
        if(program.getBlacklist() != null && !program.getBlacklist().getExclusion().isEmpty())
        {
            for (ProgramBlacklist.Exclusion exclusion : program.getBlacklist().getExclusion())
            {
                String hgvsProtein = exclusion.getHGVSP() != null ? exclusion.getHGVSP() : "";
                int minCodon = exclusion.getMinCodon() != null ? exclusion.getMinCodon().intValue() : -1;

                String chromosome = "";
                long position = 0;
                String ref = "";
                String alt = "";

                if(exclusion.getPosition() != null)
                {
                    String[] varDetails = exclusion.getPosition().split(":");

                    if(varDetails.length == 4)
                    {
                        chromosome = varDetails[0];
                        position = Long.parseLong(varDetails[1]);
                        ref = varDetails[2];
                        alt = varDetails[3];
                    }
                }

                final String gene = exclusion.getGene().getName();

                VariantFilter filter = new VariantFilter(gene, "", chromosome, position,
                        ref, alt, NONSENSE_OR_FRAMESHIFT, hgvsProtein, "", "", minCodon);

                List<VariantFilter> filters = mBlacklistFilters.get(gene);

                if(filters == null)
                {
                    filters = Lists.newArrayList();
                    mBlacklistFilters.put(gene, filters);
                }

                filters.add(filter);
            }
        }

        if(program.getWhitelist() != null && !program.getWhitelist().getVariantOrDbSNP().isEmpty())
        {
            // add to generic filter collection if to be used
        }

        return true;
    }

    public void addExternalFilters(final List<VariantFilter> allFilters)
    {
        // split into white and black list based on the coding effect
        for(VariantFilter filter : allFilters)
        {
            List<VariantFilter> filters = null;

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

    List<EligibilityReport> processVcfFile(final String sampleId, final VCFFileReader reader, boolean usesIndex)
    {
        mReports.clear();

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

        return mReports;
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
        EligibilityReport.MatchType matchType = NONE;

        // first check the transcript

        for (int i = 0; i < sampleAnnotations.size(); ++i)
        {
            final SnpEffAnnotation snpEff = sampleAnnotations.get(i);

            if(!snpEff.allele().equals(refGenotype.getAllele(1).getBaseString()))
                continue; // skip alleles for the second sample (ie the tumor)

            if (!snpEff.isTranscriptFeature())
                continue;

            if (region != null)
            {
                if(!region.transcriptID().equals(snpEff.transcript()))
                    continue;
            }
            else
            {
                if(!mPanelTranscripts.contains(snpEff.transcript()))
                    continue;
            }

            final String gene = snpEff.gene();
            CodingEffect codingEffect = CodingEffect.effect(gene, snpEff.consequences());

            final String varId = variant.getID();
            final String transcriptId = snpEff.transcript();
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
                    LOGGER.debug("match found: gene({} {}) var({}:{}) ref({}) alt({}) on effect({})",
                            gene, transcriptId, chromosome, position, ref, alt, effects);

                    matchType = REQUIRED_EFFECT;
                    break;
                }
            }

            boolean matchesFilter = false;

            if (matchType == REQUIRED_EFFECT && !mBlacklistFilters.isEmpty())
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
                        matchesFilter = filter.blacklistMatch(gene, chromosome, position, ref, alt, proteinPosition);

                        if(matchesFilter)
                        {
                            if(isBenign(filter.ClinvarSignificance) || filter.ClinvarSignificance.isEmpty())
                            {
                                LOGGER.debug("gene({}) var({}:{}:{}) ref({}) alt({}) protein({}) blacklisted",
                                        gene, varId, chromosome, position, ref, alt, hgvsProtein);

                                return;
                            }
                        }
                    }
                }
            }

            if (matchType == NONE && !mWhitelistFilters.isEmpty())
            {
                List<VariantFilter> filters = mWhitelistFilters.get(gene);

                if(filters != null)
                {
                    for (final VariantFilter filter : filters)
                    {
                        if(filter.whitelistMatch(gene, chromosome, position, ref, alt, codingEffect, hgvsProtein))
                        {
                            LOGGER.debug("match found: gene({} {}) var({}:{}:{}) ref({}) alt({}) hgvsProtein({}) whitelisted",
                                    gene, transcriptId, varId, chromosome, position, ref, alt, hgvsProtein);
                            matchType = WHITELIST;
                            matchesFilter = true;
                            break;
                        }
                    }
                }
            }

            if(matchType == NONE)
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

            EligibilityReport report = ImmutableEligibilityReport.builder()
                    .sampleId(sampleId)
                    .program(mName)
                    .matchType(matchType)
                    .id(variant.getID())
                    .genes(gene)
                    .transcriptId(transcriptId)
                    .chrom(chromosome)
                    .pos(position)
                    .ref(ref)
                    .alts(alt)
                    .effects(effects)
                    .codingEffect(codingEffect)
                    .annotations(annotationsStr)
                    .hgvsProtein(hgvsProtein)
                    .hgvsCoding(hgvsCoding)
                    .isHomozygous(isHomozygous)
                    .phredScore(phredScore)
                    .hasDepthInfo(hasDepthInfo)
                    .germlineAltCount(germlineAltCount)
                    .germlineReadDepth(germlineReadDepth)
                    .tumorAltCount(tumorAltCount)
                    .tumorReadDepth(tumorReadDepth)
                    .condonInfo(codonInfo)
                    .matchesClinvarFilter(matchesFilter)
                    .build();

            mReports.add(report);
        }
    }

    private static List<Integer> proteinPosition(final SnpEffAnnotation annotation)
    {
        return Arrays.stream(annotation.aaPosAndLength().split("/"))
                .filter(s -> !s.isEmpty())
                .map(Integer::parseInt)
                .collect(Collectors.toList());
    }

    public static boolean matchesWhitelistGeneProtein(ProgramWhitelist.Variant geneProtein,
            final VariantContext context, final SnpEffAnnotation annotation)
    {
        if(!geneProtein.getGene().getName().equals(annotation.gene()))
            return false;

        if (geneProtein.getHGVSP() != null && !annotation.hgvsProtein().isEmpty()
                && geneProtein.getHGVSP().equals(annotation.hgvsProtein().replaceFirst("^p\\.", "")))
        {
            LOGGER.debug("variant({}) found in blacklist HGVSP({})", context.getID(), geneProtein.getHGVSP());
            return true;
        }

        return false;
    }

    public static boolean matchesBlacklistExclusion(final ProgramBlacklist.Exclusion blacklist, final VariantContext context,
            final SnpEffAnnotation annotation) {

        if (blacklist.getHGVSP() != null && !annotation.hgvsProtein().isEmpty()
        && blacklist.getHGVSP().equals(annotation.hgvsProtein().replaceFirst("^p\\.", "")))
        {
            LOGGER.debug("variant({}) found in blacklist HGVSP({})", context.getID(), blacklist.getHGVSP());
            return true;
        }

        if (blacklist.getHGVSC() != null && !annotation.hgvsCoding().isEmpty()
        && blacklist.getHGVSC().equals(annotation.hgvsCoding().replaceFirst("^c\\.", ""))) {

            LOGGER.debug("variant({}) found in blacklist HGVSC({})", context.getID(), blacklist.getHGVSC());

            return true;
        }

        final List<Integer> proteinPositions = proteinPosition(annotation);
        if (blacklist.getMinCodon() != null && !proteinPositions.isEmpty()
                && blacklist.getMinCodon().intValue() <= proteinPositions.get(0))
        {

            LOGGER.debug("variant({}) found in blacklist minCodon({})", context.getID(), blacklist.getMinCodon());
            return true;
        }

        if(blacklist.getPosition() != null && atPosition(context, blacklist.getPosition()))
        {
            LOGGER.debug("variant({}) found in blacklist postition({})", context.getID(), blacklist.getPosition());
            return true;
        }

        return false;
    }

    private static boolean atPosition(final VariantContext v, final String position)
    {
        return position.equals(v.getContig() + ":" + v.getStart());
    }


    private void initialiseGeneData()
    {
        mGenesByChromosomeMap = HmfGenePanelSupplier.allGenesPerChromosomeMap37();

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

    /*  OLD BLACKLIST PREDICATE

        public boolean test(final VariantModel variantModel)
        {
            for (final SnpEffAnnotation annotation : variantModel.sampleAnnotations())
            {
                final boolean transcriptMatches = transcripts.contains(annotation.transcript());
                if (transcriptMatches)
                {
                    for (ProgramBlacklist.Exclusion exclusion : blacklist)
                    {
                        if (matchesBlacklistExclusion(exclusion, variantModel.context(), annotation))
                        {
                            return true;
                        }
                    }
                }
            }
            return false;
        }


        public static List<Integer> proteinPosition(@NotNull final SnpEffAnnotation annotation)
        {
            return Arrays.stream(annotation.aaPosAndLength().split("/"))
                    .filter(s -> !s.isEmpty())
                    .map(Integer::parseInt)
                    .collect(Collectors.toList());
        }

    */

    /* OLD WHITELIST PREDICATE

        public boolean test(final VariantModel variantModel)
        {
            if(mWhitelist == null)
                return false;

            for (final SnpEffAnnotation annotation : variantModel.sampleAnnotations())
            {
                for (final Object variantOrDbSNP : mWhitelist.getVariantOrDbSNP())
                {
                    if (variantOrDbSNP instanceof ProgramWhitelist.Variant)
                    {
                        final ProgramWhitelist.Variant whitelistVar = (ProgramWhitelist.Variant) variantOrDbSNP;

                        if(matchesWhitelistGeneProtein(whitelistVar, variantModel.context(), annotation))
                        {
                            return true;
                        }
                    }
                    else if (variantOrDbSNP instanceof String)
                    {
                        if(matchesWhitelistDbSNPId((String) variantOrDbSNP, variantModel.context()))
                        {
                            return true;
                        }
                    }
                }
            }

            return false;
        }

        public static boolean matchesWhitelistDbSNPId(final String dnSNPId, final VariantContext variant)
        {
            Set<String> varDbSNPList = Lists.newArrayList(variant.getID()
                    .split(","))
                    .stream().filter(s -> s.startsWith("rs"))
                    .collect(Collectors.toSet());

            if(varDbSNPList.contains(dnSNPId))
            {
                LOGGER.debug("variant({}) matched in whitelist rs DB Ids list({})", variant.getID());
                return true;
            }
            else
            {
                return false;
            }
        }
    */

}
