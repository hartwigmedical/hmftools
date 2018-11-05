package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.GENE_TRANSCRIPT;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.HOTSPOT_LOCATION;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.NONE;
import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.WHITELIST;
import static com.hartwig.hmftools.bachelor.EligibilityReport.ReportType.GERMLINE_DELETION;
import static com.hartwig.hmftools.bachelor.EligibilityReport.ReportType.SOMATIC_DELETION;
import static com.hartwig.hmftools.bachelor.EligibilityReport.ReportType.SOMATIC_DISRUPTION;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.bachelor.predicates.BlacklistPredicate;
import com.hartwig.hmftools.bachelor.predicates.WhitelistPredicate;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.region.HmfExonRegion;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import nl.hartwigmedicalfoundation.bachelor.GeneIdentifier;
import nl.hartwigmedicalfoundation.bachelor.HotspotLocation;
import nl.hartwigmedicalfoundation.bachelor.OtherEffect;
import nl.hartwigmedicalfoundation.bachelor.Program;
import nl.hartwigmedicalfoundation.bachelor.ProgramPanel;
import nl.hartwigmedicalfoundation.bachelor.SnpEffect;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

class BachelorEligibility {

    private static final Logger LOGGER = LogManager.getLogger(BachelorEligibility.class);
    private static final SortedSetMultimap<String, HmfTranscriptRegion> ALL_GENES_BY_CHROMOSOME = HmfGenePanelSupplier.allGenesPerChromosomeMap37();
    private static final Map<String, HmfTranscriptRegion> ALL_GENES = makeGeneNameMap();
    private static final Map<String, HmfTranscriptRegion> ALL_TRANSCRIPT_IDS = makeTranscriptMap();

    private final List<BachelorProgram> programs = Lists.newArrayList();
    private final Set<HmfTranscriptRegion> variantLocationsToQuery = Sets.newHashSet();

    @NotNull
    private static Map<String, HmfTranscriptRegion> makeGeneNameMap()
    {
        final Map<String, HmfTranscriptRegion> result = Maps.newHashMap();

        for (final HmfTranscriptRegion region : ALL_GENES_BY_CHROMOSOME.values())
        {
            result.put(region.gene(), region);
        }
        return result;
    }

    @NotNull
    private static Map<String, HmfTranscriptRegion> makeTranscriptMap()
    {
        final Map<String, HmfTranscriptRegion> result = Maps.newHashMap();

        for (final HmfTranscriptRegion region : ALL_GENES_BY_CHROMOSOME.values())
        {
            result.put(region.transcriptID(), region);
        }
        return result;
    }

    private BachelorEligibility() {
    }

    @NotNull
    static BachelorEligibility fromMap(@NotNull Map<String, Program> input)
    {
        final BachelorEligibility result = new BachelorEligibility();

        for (final Program program : input.values())
        {
            final Multimap<String, String> geneToEnsemblMap = HashMultimap.create();

            program.getPanel()
                    .stream()
                    .map(ProgramPanel::getGene)
                    .flatMap(Collection::stream)
                    .forEach(g -> geneToEnsemblMap.put(g.getName(), g.getEnsembl()));

            // process variants from vcf
            final List<Predicate<VariantModel>> panelPredicates = Lists.newArrayList();

            List<String> requiredEffects = Lists.newArrayList();
            List<String> panelTranscripts = Lists.newArrayList();
            List<HotspotLocation> hotspots = Lists.newArrayList();

            for (final ProgramPanel panel : program.getPanel())
            {
                final List<GeneIdentifier> genes = panel.getGene();

                // take up a collection of the effects to search for
                hotspots = panel.getHotspot();
                requiredEffects = panel.getSnpEffect().stream().map(SnpEffect::value).collect(Collectors.toList());
                panelTranscripts = genes.stream().map(GeneIdentifier::getEnsembl).collect(Collectors.toList());


                /*
                final List<String> effects = requiredEffects;

                final Predicate<VariantModel> panelPredicate = v -> genes.stream()
                        .anyMatch(p -> v.sampleAnnotations()
                                .stream()
                                .anyMatch(a -> a.transcript().equals(p.getEnsembl()) && effects.stream()
                                        .anyMatch(x -> a.effects().contains(x))));

                final Predicate<VariantModel> effectsPredicate = v -> effects.stream()
                        .anyMatch(p -> v.sampleAnnotations()
                                .stream()
                                .anyMatch(a -> a.effects().contains(p)));

                Predicate<VariantModel> transcriptPredicate = v -> genes.stream()
                        .anyMatch(p -> v.sampleAnnotations()
                                .stream()
                                .anyMatch(a -> a.transcript().equals(p.getEnsembl())));

                transcriptPredicate = transcriptPredicate.and(effectsPredicate);

                Predicate<VariantModel> hotspotPredicate = v -> hotspots.stream()
                        .anyMatch(p -> v.context().getContig().equals(p.getChromosome())
                                && v.context().getStart() == p.getPosition().intValue()
                                && v.context().getReference().getBaseString().equals(p.getRef())
                                && v.context().getAlleles().size() >= 2 && v.context().getAlleles().get(1).getBaseString().equals(p.getAlt()));

                hotspotPredicate = hotspotPredicate.and(effectsPredicate);

                final Predicate<VariantModel> combinedPredicate = transcriptPredicate.or(hotspotPredicate);

                panelPredicates.add(combinedPredicate);

                panelPredicates.add(effectsPredicate);
                */

                // update query targets
                for (final GeneIdentifier g : genes)
                {
                    final HmfTranscriptRegion region = ALL_TRANSCRIPT_IDS.get(g.getEnsembl());

                    if (region == null)
                    {
                        final HmfTranscriptRegion namedRegion = ALL_GENES.get(g.getName());

                        if (namedRegion == null)
                        {
                            LOGGER.warn("Program {} gene {} non-canonical transcript {} couldn't find region, transcript will be skipped",
                                    program.getName(), g.getName(), g.getEnsembl());

                            // just skip this gene for now
                        }
                        else
                        {
                            result.variantLocationsToQuery.add(namedRegion);
                        }
                    }
                    else
                    {
                        result.variantLocationsToQuery.add(region);
                    }
                }
            }

            // final Predicate<VariantModel> inPanel = v -> panelPredicates.stream().anyMatch(p -> p.test(v));
            final Predicate<VariantModel> inPanel = v -> true; // manually checked for each variant since too difficult to express as a predicate

            final Predicate<VariantModel> inBlacklist = new BlacklistPredicate(geneToEnsemblMap.values(), program.getBlacklist());
            final Predicate<VariantModel> inWhitelist = new WhitelistPredicate(geneToEnsemblMap, program.getWhitelist());
            // final Predicate<VariantModel> snvPredicate = v -> inPanel.test(v) ? !inBlacklist.test(v) : inWhitelist.test(v);
            final Predicate<VariantModel> snvPredicate = v -> !inBlacklist.test(v);

            BachelorProgram bachelorProgram = new BachelorProgram(program.getName(), snvPredicate, inWhitelist, requiredEffects, panelTranscripts, hotspots);

            result.programs.add(bachelorProgram);
        }

        return result;
    }

    @NotNull
    private Collection<EligibilityReport> processVariant(final VariantContext variant, final String patient, final String sample, final EligibilityReport.ReportType type)
    {
        if (variant.isFiltered())
            return Collections.emptyList();

        // we will skip when an ALT is not present in the sample
        final Genotype genotype = variant.getGenotype(sample);

        if (genotype == null || !(genotype.isHomVar() || genotype.isHet()))
        {
            return Collections.emptyList();
        }

        // gather up the relevant alleles
        VariantModel sampleVariant = new VariantModel(sample, variant);

        // apply the all relevant tests to see if this program has been matched
        final List<String> matchingPrograms = programs.stream()
                .filter(program -> program.vcfProcessor().test(sampleVariant))
                .map(BachelorProgram::name)
                .collect(Collectors.toList());

        List<EligibilityReport> reportList = Lists.newArrayList();

        // search the list of annotations for the correct allele and transcript ID to write to the result file
        for (BachelorProgram program : programs)
        {
            final List<String> requiredEffects = program.requiredEffects();
            final List<String> geneTranscripts = program.panelTranscripts();
            final List<HotspotLocation> hotspots = program.hotspots();

            // if (!program.vcfProcessor().test(sampleVariant))
            //     continue;

            // check the sub-conditions now - hotspot locations and gene-transcript IDs

            EligibilityReport.MatchType matchType = NONE;

            // first check the transcript
            SnpEffAnnotation relevantSnpEff = null;
            String annotationsStr = "";

            for (int i = 0; i < sampleVariant.sampleAnnotations().size(); ++i)
            {
                final SnpEffAnnotation snpEff = sampleVariant.sampleAnnotations().get(i);

                if (!snpEff.isTranscriptFeature())
                    continue;

                if (geneTranscripts.contains(snpEff.transcript()))
                {
                    for (String requiredEffect : requiredEffects)
                    {
                        if (snpEff.effects().contains(requiredEffect))
                        {
                            LOGGER.debug("match found: program({}): var({}:{}) ref({}) alt({}) on effect({}) and transcript({})",
                                    program.name(), sampleVariant.context().getContig(), sampleVariant.context().getStart(),
                                    sampleVariant.context().getReference().getBaseString(), sampleVariant.context().getAlleles().get(1).getBaseString(),
                                    snpEff.effects(), snpEff.transcript());

                            matchType = GENE_TRANSCRIPT;
                            relevantSnpEff = snpEff;
                            annotationsStr = sampleVariant.rawAnnotations().get(i);
                            break;
                        }
                    }
                }

                if (matchType == GENE_TRANSCRIPT)
                    break;
            }

            // then check the hotspot location

            if(matchType == NONE && program.whitelist().test(sampleVariant))
            {
                matchType = WHITELIST;
            }

            if(matchType == NONE)
            {
                for (final HotspotLocation hotspot : hotspots)
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
                            program.name(), sampleVariant.context().getContig(), sampleVariant.context().getStart(),
                            sampleVariant.context().getReference().getBaseString(), sampleVariant.context() .getAlleles() .get(1) .getBaseString());
                }
            }

            if(matchType == HOTSPOT_LOCATION || matchType == WHITELIST)
            {
                // select the first relevant feature
                for (int i = 0; i < sampleVariant.sampleAnnotations().size(); ++i)
                {
                    final SnpEffAnnotation snpEff = sampleVariant.sampleAnnotations().get(i);

                    if (!snpEff.isTranscriptFeature())
                        continue;

                    relevantSnpEff = snpEff;
                    annotationsStr = sampleVariant.rawAnnotations().get(i);
                    break;
                }
            }

            if (matchType == NONE)
            {
                continue;
            }

            boolean isHomozygous = variant.getGenotype(0).isHom();
            int phredScore = variant.getGenotype(0).getPL().length >= 1 ? variant.getGenotype(0).getPL()[0] : 0;

            EligibilityReport report = ImmutableEligibilityReport.builder()
                    .patient(patient)
                    .source(type)
                    .program(program.name())
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
                    .build();

            reportList.add(report);
        }

        /*
        if (!reportList.isEmpty()) {
            LOGGER.debug("writing {} matched reports", reportList.size());
        }
        */

        return reportList;
    }

    @NotNull
    Collection<EligibilityReport> processVCF(final String patient, final String sample, final EligibilityReport.ReportType type, final VCFFileReader reader)
    {
        final List<EligibilityReport> results = Lists.newArrayList();

        for (final HmfTranscriptRegion region : variantLocationsToQuery)
        {

            // LOGGER.debug("chromosome({} start={} end={})", region.chromosome(), (int) region.geneStart(), (int) region.geneEnd());

            final CloseableIterator<VariantContext> query =
                    reader.query(region.chromosome(), (int) region.geneStart(), (int) region.geneEnd());

            while (query.hasNext()) {
                final VariantContext variant = query.next();
                // LOGGER.debug("patient({}) sample({}) processing variant({})", patient, sample, variant.getID());
                results.addAll(processVariant(variant, patient, sample, type));
            }
            query.close();
        }

        return results;
    }

}
