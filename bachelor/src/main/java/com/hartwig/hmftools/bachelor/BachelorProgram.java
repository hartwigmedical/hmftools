package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.GENE_TRANSCRIPT;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.HOTSPOT_LOCATION;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.NONE;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.WHITELIST;

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
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import nl.hartwigmedicalfoundation.bachelor.GeneIdentifier;
import nl.hartwigmedicalfoundation.bachelor.HotspotLocation;
import nl.hartwigmedicalfoundation.bachelor.Program;
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
    private List<HotspotLocation> mHotspots;
    private Predicate<VariantModel> mWhiteList;
    private Predicate<VariantModel> mBlackList;
    private boolean mHasWhitelist;
    private boolean mHasBlacklist;

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

        mBlackList = null;
        mWhiteList = null;
        mHasWhitelist = false;
        mHasBlacklist = false;

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

        if(program.getBlacklist() != null && !program.getBlacklist().getExclusion().isEmpty())
        {
            final Predicate<VariantModel> inBlacklist = new BlacklistPredicate(geneToEnsemblMap.values(), program.getBlacklist());
            mBlackList = v -> !inBlacklist.test(v);
            mHasBlacklist = true;
        }

        if(program.getWhitelist() != null && !program.getWhitelist().getVariantOrDbSNP().isEmpty())
        {
            mWhiteList = new WhitelistPredicate(geneToEnsemblMap, program.getWhitelist());
            mHasWhitelist = program.getWhitelist() != null && !program.getWhitelist().getVariantOrDbSNP().isEmpty();
        }

        return true;
    }

    public String name() { return mName; }

    List<EligibilityReport> processVCF(final String sampleId, final VCFFileReader reader)
    {
        mReports.clear();

        for (final HmfTranscriptRegion region : mTranscriptRegions)
        {
            // LOGGER.debug("chromosome({} start={} end={})", region.chromosome(), (int) region.geneStart(), (int) region.geneEnd());

            final CloseableIterator<VariantContext> query = reader.query(region.chromosome(), (int) region.geneStart(), (int) region.geneEnd());

            while (query.hasNext())
            {
                final VariantContext variant = query.next();
                processVariant(variant, sampleId, region);
            }

            query.close();
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
        SnpEffAnnotation relevantSnpEff = null;
        String annotationsStr = "";

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

            for (String requiredEffect : mRequiredEffects)
            {
                if (snpEff.effects().contains(requiredEffect))
                {
                    LOGGER.debug("match found: program({}): var({}:{}) ref({}) alt({}) on effect({}) and transcript({})",
                            mName, variant.getContig(), variant.getStart(),
                            variant.getReference().getBaseString(), variant.getAlleles().get(1).getBaseString(),
                            snpEff.effects(), snpEff.transcript());

                    matchType = GENE_TRANSCRIPT;
                    relevantSnpEff = snpEff;
                    annotationsStr = SnpEffAnnotationFactory.rawAnnotations(variant).get(i);
                    break;
                }
            }

            if (matchType == GENE_TRANSCRIPT)
                break;
        }

        if(matchType == NONE && mHasWhitelist)
        {
            VariantModel sampleVariant = new VariantModel(refGenotype.getSampleName(), variant);

            if(mWhiteList.test(sampleVariant))
                matchType = WHITELIST;
        }

        // then check the hotspot location
        if(matchType == NONE)
        {
            for (final HotspotLocation hotspot : mHotspots)
            {
                if (variant.getStart() != hotspot.getPosition().intValue() || !variant.getContig().equals(hotspot.getChromosome()))
                    continue;

                if (!variant.getReference().getBaseString().equals(hotspot.getRef())
                        || variant.getAlleles().size() < 2 || !variant.getAlleles().get(1).getBaseString().equals(hotspot.getAlt()))
                {
                    continue;
                }

                matchType = HOTSPOT_LOCATION;

                LOGGER.debug("match found: program({}): var({}:{}) ref({}) alt({}) on hotspot location",
                        mName, variant.getContig(), variant.getStart(),
                        variant.getReference().getBaseString(), variant .getAlleles() .get(1) .getBaseString());
            }
        }

        if(matchType == HOTSPOT_LOCATION || matchType == WHITELIST)
        {
            // select the first relevant feature
            for (int i = 0; i < sampleAnnotations.size(); ++i)
            {
                final SnpEffAnnotation snpEff = sampleAnnotations.get(i);

                if (!snpEff.isTranscriptFeature())
                    continue;

                relevantSnpEff = snpEff;
                annotationsStr = SnpEffAnnotationFactory.rawAnnotations(variant).get(i);
                break;
            }
        }

        // check blacklistings
        if(mHasBlacklist)
        {
            VariantModel sampleVariant = new VariantModel(refGenotype.getSampleName(), variant);

            if(!mBlackList.test(sampleVariant))
            {
                LOGGER.debug("var({}:{}) ref({}) alt({}) blacklisted",
                        variant.getContig(), variant.getStart(), variant.getReference().getBaseString(),
                        variant .getAlleles().get(1).getBaseString());
                return;
            }
        }

        if (matchType == NONE)
            return;

        boolean isHomozygous = refGenotype.isHom();
        int phredScore = refGenotype.getPL().length >= 1 ? refGenotype.getPL()[0] : 0;

        int germlineAltCount = refGenotype.getAD()[1];
        int germlineReadDepth = refGenotype.getDP();

        int tumorAltCount = 0;
        int tumorReadDepth = 0;

        if(variant.getGenotypes().size() >= 2)
        {
            final Genotype tumorGenotype = variant.getGenotype(1);
            int[] alleleData = tumorGenotype.getAD();
            tumorAltCount = alleleData[1];
            tumorReadDepth = tumorGenotype.getDP();
        }

        final String codonInfo = relevantSnpEff.aaPosAndLength();

        EligibilityReport report = ImmutableEligibilityReport.builder()
                .sampleId(sampleId)
                .source(EligibilityReport.ReportType.GERMLINE_MUTATION)
                .program(mName)
                .matchType(matchType)
                .id(variant.getID())
                .genes(relevantSnpEff.gene())
                .transcriptId(relevantSnpEff.transcript())
                .chrom(variant.getContig())
                .pos(variant.getStart())
                .ref(variant.getReference().toString())
                .alts(relevantSnpEff.allele())
                .effects(relevantSnpEff.effects())
                .annotations(annotationsStr)
                .hgvsProtein(relevantSnpEff.hgvsProtein())
                .hgvsCoding(relevantSnpEff.hgvsCoding())
                .isHomozygous(isHomozygous)
                .phredScore(phredScore)
                .germlineAltCount(germlineAltCount)
                .germlineReadDepth(germlineReadDepth)
                .tumorAltCount(tumorAltCount)
                .tumorReadDepth(tumorReadDepth)
                .condonInfo(codonInfo)
                .build();

        mReports.add(report);
    }

    public boolean hasWhiteList() { return mHasWhitelist; }
    public Predicate<VariantModel> whitelist() { return mWhiteList; }

    public List<String> requiredEffects() { return mRequiredEffects; }
    public List<String> panelTranscripts() { return mPanelTranscripts; }
    public List<HotspotLocation> hotspots() { return mHotspots; }

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
