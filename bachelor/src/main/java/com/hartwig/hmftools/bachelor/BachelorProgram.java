package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.GENE_TRANSCRIPT;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.HOTSPOT_LOCATION;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.NONE;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.WHITELIST;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.bachelor.predicates.BlacklistPredicate;
import com.hartwig.hmftools.bachelor.predicates.WhitelistPredicate;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import nl.hartwigmedicalfoundation.bachelor.GeneIdentifier;
import nl.hartwigmedicalfoundation.bachelor.HotspotLocation;
import nl.hartwigmedicalfoundation.bachelor.Program;
import nl.hartwigmedicalfoundation.bachelor.ProgramBlacklist;
import nl.hartwigmedicalfoundation.bachelor.ProgramPanel;
import nl.hartwigmedicalfoundation.bachelor.SnpEffect;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class BachelorProgram
{
    private String mName;

    private List<String> mRequiredEffects;
    private List<String> mPanelTranscripts;
    private List<HotspotLocation> mHotspots;

    private List<VariantFilter> mWhitelistFilters;
    private List<VariantFilter> mBlacklistFilters;

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
        mHotspots = null;

        mWhitelistFilters = Lists.newArrayList();
        mBlacklistFilters = Lists.newArrayList();

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
        mHotspots = panel.getHotspot();
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

                if(exclusion.getPosition() != null)
                {
                    String[] chrPos = exclusion.getPosition().split(":");

                    if(chrPos.length == 2)
                    {
                        chromosome = chrPos[0];
                        position = Long.parseLong(chrPos[1]);
                    }
                }

                VariantFilter filter = new VariantFilter(exclusion.getGene().getName(), "", chromosome, position,
                        "", "", NONSENSE_OR_FRAMESHIFT, hgvsProtein, "", minCodon);

                mBlacklistFilters.add(filter);
            }
        }

        if(program.getWhitelist() != null && !program.getWhitelist().getVariantOrDbSNP().isEmpty())
        {
            // add to generic filter collection if to be used
        }

        return true;
    }

    public void addExternalFilters(final List<VariantFilter> filters)
    {
        // split into white and black list based on the coding effect
        mBlacklistFilters.addAll(filters.stream()
                .filter(x -> x.Effect == NONSENSE_OR_FRAMESHIFT || x.Effect == CodingEffect.SPLICE)
                .collect(Collectors.toList()));

        mWhitelistFilters.addAll(filters.stream()
                .filter(x -> x.Effect != NONSENSE_OR_FRAMESHIFT && x.Effect != CodingEffect.SPLICE)
                .collect(Collectors.toList()));
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
            final String ref = variant.getReference().toString();
            final String alt = snpEff.allele();
            final String effects = snpEff.effects();
            final String hgvsProtein = snpEff.hgvsProtein();
            final String hgvsCoding = snpEff.hgvsCoding();

            for (String requiredEffect : mRequiredEffects)
            {
                if (effects.contains(requiredEffect))
                {
                    LOGGER.debug("match found: gene({} {}) var({}:{}) ref({}) alt({}) on effect({})",
                            gene, transcriptId, chromosome, position, ref, alt, effects);

                    matchType = GENE_TRANSCRIPT;
                    break;
                }
            }

            if (matchType == GENE_TRANSCRIPT && !mBlacklistFilters.isEmpty())
            {
                for(final VariantFilter filter : mBlacklistFilters)
                {
                    if(!filter.Gene.equals(gene))
                        continue;

                    if(!filter.HgvsProteinCodon.isEmpty())
                    {
                        if(filter.HgvsProteinCodon.equals(hgvsProtein))
                        {
                            LOGGER.debug("gene({} {}) var({}:{}) ref({}) alt({}) blacklisted on hgvsProtein({})",
                                    gene, transcriptId, chromosome, position, ref, alt, hgvsProtein);
                            return;
                        }

                        continue;
                    }

                    if(!filter.DBSnpId.isEmpty())
                    {
                        if(varId.contains(filter.DBSnpId))
                        {
                            LOGGER.debug("gene({} {}) var({}:{}) ref({}) alt({}) blacklisted on DBSnpId({})",
                                    gene, transcriptId, chromosome, position, ref, alt, varId);
                            return;
                        }

                        continue;
                    }

                    if(filter.MinCodon >= 0)
                    {
                        final List<Integer> proteinPositions = proteinPosition(snpEff);

                        if(!proteinPositions.isEmpty() && filter.MinCodon <= proteinPositions.get(0))
                        {
                            LOGGER.debug("gene({} {}) var({}:{}) ref({}) alt({}) blacklisted on minCodon({})",
                                    gene, transcriptId, chromosome, position, ref, alt, proteinPositions.get(0));
                            return;
                        }
                    }

                    if(filter.Position == position && filter.Ref.equals(ref) && filter.Alt.equals(alt))
                    {
                        LOGGER.debug("gene({} {}) var({}:{}) ref({}) alt({}) blacklisted on position and ref/alt",
                                gene, transcriptId, chromosome, position, ref, alt);
                        return;
                    }

                }
            }

            if (matchType == NONE && !mWhitelistFilters.isEmpty())
            {
                for(final VariantFilter filter : mWhitelistFilters)
                {
                    if(!filter.Gene.equals(gene))
                        continue;

                    if(!filter.HgvsProteinCodon.isEmpty() && filter.HgvsProteinCodon.equals(hgvsProtein))
                    {
                        LOGGER.debug("match found: gene({} {}) var({}:{}) ref({}) alt({}) on hgvsProtein({}) whitelist",
                                gene, transcriptId, chromosome, position, ref, alt, hgvsProtein);
                        matchType = WHITELIST;
                        break;
                    }

                    if(!filter.DBSnpId.isEmpty() && varId.contains(filter.DBSnpId))
                    {
                        LOGGER.debug("match found: gene({} {}) var({}:{}) ref({}) alt({}) on DBSnpId({}) whitelist",
                                gene, transcriptId, chromosome, position, ref, alt, varId);
                        matchType = WHITELIST;
                        break;
                    }
                }
            }

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

            if(matchType == NONE)
                return;

            String annotationsStr = SnpEffAnnotationFactory.rawAnnotations(variant).get(i);

            boolean isHomozygous = refGenotype.isHom();
            int phredScore = refGenotype.getPL().length >= 1 ? refGenotype.getPL()[0] : 0;

            int germlineAltCount = refGenotype.getAD()[1];
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

}
