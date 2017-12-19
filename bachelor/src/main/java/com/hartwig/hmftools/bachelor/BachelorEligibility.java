package com.hartwig.hmftools.bachelor;

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
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.region.hmfslicer.HmfExonRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.hmfslicer.HmfGenePanelSupplier;

import nl.hartwigmedicalfoundation.bachelor.GeneIdentifier;
import nl.hartwigmedicalfoundation.bachelor.OtherEffect;
import nl.hartwigmedicalfoundation.bachelor.Program;
import nl.hartwigmedicalfoundation.bachelor.ProgramBlacklist;
import nl.hartwigmedicalfoundation.bachelor.ProgramPanel;
import nl.hartwigmedicalfoundation.bachelor.ProgramWhitelist;
import nl.hartwigmedicalfoundation.bachelor.SnpEffect;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

class BachelorEligibility {

    private static final double MAX_COPY_NUMBER_FOR_LOSS = 0.5;

    private static final Logger LOGGER = LogManager.getLogger(BachelorEligibility.class);
    private static final SortedSetMultimap<String, HmfGenomeRegion> allGenesByChromosomeMap = HmfGenePanelSupplier.allGeneMap();
    private static final Map<String, HmfGenomeRegion> allGenesMap = makeGeneNameMap();
    private static final Map<String, HmfGenomeRegion> allTranscriptsMap = makeTranscriptMap();

    private final Map<String, BachelorProgram> programs = Maps.newHashMap();
    private final Map<String, HmfGenomeRegion> querySet = Maps.newHashMap();

    private static Map<String, HmfGenomeRegion> makeGeneNameMap() {
        final Map<String, HmfGenomeRegion> result = Maps.newHashMap();
        for (final HmfGenomeRegion region : allGenesByChromosomeMap.values()) {
            result.put(region.gene(), region);
        }
        return result;
    }

    private static Map<String, HmfGenomeRegion> makeTranscriptMap() {
        final Map<String, HmfGenomeRegion> result = Maps.newHashMap();
        for (final HmfGenomeRegion region : allGenesByChromosomeMap.values()) {
            result.put(region.transcriptID(), region);
        }
        return result;
    }

    private BachelorEligibility() {
    }

    private static boolean atPosition(final VariantContext v, final String position) {
        // TODO: robust enough check?
        return position.equals(v.getContig() + ":" + v.getStart());
    }

    static BachelorEligibility fromMap(final Map<String, Program> input) {
        final BachelorEligibility result = new BachelorEligibility();
        boolean iterateAllGenes = false;

        for (final Program program : input.values()) {

            final Multimap<String, String> geneToEnsemblMap = HashMultimap.create();
            program.getPanel()
                    .stream()
                    .map(ProgramPanel::getGene)
                    .flatMap(Collection::stream)
                    .forEach(g -> geneToEnsemblMap.put(g.getName(), g.getEnsembl()));

            // process copy number sections
            final List<Predicate<GeneCopyNumber>> cnvPredicates = Lists.newArrayList();
            for (final ProgramPanel panel : program.getPanel()) {
                final boolean allGene = panel.getAllGenes() != null;
                final List<GeneIdentifier> genes = panel.getGene();

                if (panel.getEffect().contains(OtherEffect.HOMOZYGOUS_DELETION)) {
                    final Predicate<GeneCopyNumber> geneCopyNumberPredicate =
                            cnv -> allGene || genes.stream().anyMatch(g -> g.getEnsembl().equals(cnv.transcriptID()));
                    cnvPredicates.add(geneCopyNumberPredicate);
                }
            }

            // process structural variant disruptions
            final List<Predicate<HmfGenomeRegion>> disruptionPredicates = Lists.newArrayList();
            for (final ProgramPanel panel : program.getPanel()) {
                final boolean allGene = panel.getAllGenes() != null;
                final List<GeneIdentifier> genes = panel.getGene();

                if (panel.getEffect().contains(OtherEffect.GENE_DISRUPTION)) {
                    final Predicate<HmfGenomeRegion> disruptionPredicate =
                            sv -> allGene || genes.stream().anyMatch(g -> g.getEnsembl().equals(sv.transcriptID()));
                    // TODO: we are matching on transcript ID here but we only have canonical transcripts in our panel file
                    disruptionPredicates.add(disruptionPredicate);
                }
            }

            // process variants from vcf
            final List<Predicate<VariantModel>> panelPredicates = Lists.newArrayList();
            for (final ProgramPanel panel : program.getPanel()) {
                final boolean allGene = panel.getAllGenes() != null;
                final List<GeneIdentifier> genes = panel.getGene();
                final List<String> effects = panel.getSnpEffect().stream().map(SnpEffect::value).collect(Collectors.toList());

                final Predicate<VariantModel> panelPredicate = v -> allGene
                        ? v.Annotations.stream().anyMatch(a -> a.Effects.stream().anyMatch(effects::contains))
                        : genes.stream()
                                .anyMatch(p -> v.Annotations.stream()
                                        .anyMatch(a -> a.Transcript.equals(p.getEnsembl()) && a.Effects.stream()
                                                .anyMatch(effects::contains)));
                panelPredicates.add(panelPredicate);

                // update query targets
                iterateAllGenes |= allGene;
                for (final GeneIdentifier g : genes) {
                    final HmfGenomeRegion region = allTranscriptsMap.get(g.getEnsembl());
                    if (region == null) {
                        final HmfGenomeRegion namedRegion = allGenesMap.get(g.getName());
                        if (namedRegion == null) {
                            LOGGER.warn("Program {} gene {} non-canonical transcript {} couldn't find region. Performance may be degraded.",
                                    program.getName(), g.getName(), g.getEnsembl());
                            iterateAllGenes = true;
                        } else {
                            result.querySet.put(namedRegion.transcriptID(), namedRegion);
                        }
                    } else {
                        result.querySet.put(region.transcriptID(), region);
                    }
                }
            }

            final Predicate<VariantModel> inPanel = v -> panelPredicates.stream().anyMatch(p -> p.test(v));

            // blacklist

            final List<ProgramBlacklist.Exclusion> blacklist =
                    program.getBlacklist() != null ? program.getBlacklist().getExclusion() : Lists.newArrayList();

            final Predicate<VariantModel> inBlacklist = v -> blacklist.stream().anyMatch(b -> {
                for (final SnpEff annotation : v.Annotations) {
                    final boolean transcriptMatches = geneToEnsemblMap.values().contains(annotation.Transcript);
                    if (transcriptMatches) {
                        if (b.getHGVSP() != null && !annotation.HGVSp.isEmpty() && b.getHGVSP().equals(annotation.HGVSp)) {
                            return true;
                        } else if (b.getHGVSC() != null && !annotation.HGVSc.isEmpty() && b.getHGVSC().equals(annotation.HGVSc)) {
                            return true;
                        } else if (b.getMinCodon() != null && !annotation.ProteinPosition.isEmpty() // TODO: stronger check here?
                                && b.getMinCodon().intValue() <= annotation.ProteinPosition.get(0)) {
                            return true;
                        } else if (b.getPosition() != null && atPosition(v.Context, b.getPosition())) {
                            return true;
                        }
                    }
                }
                return false;
            });

            // whitelist

            final Multimap<String, String> whitelist = HashMultimap.create();
            final Set<String> dbSNP = Sets.newHashSet();
            if (program.getWhitelist() != null) {
                for (final Object o : program.getWhitelist().getVariantOrDbSNP()) {
                    if (o instanceof ProgramWhitelist.Variant) {
                        final ProgramWhitelist.Variant v = (ProgramWhitelist.Variant) o;
                        for (final String transcript : geneToEnsemblMap.get(v.getGene().getName())) {
                            whitelist.put(transcript, v.getHGVSP());
                        }
                    } else if (o instanceof String) {
                        dbSNP.add((String) o);
                    }
                }
            }
            final Predicate<VariantModel> inWhitelist = v -> v.dbSNP.stream().anyMatch(dbSNP::contains) || v.Annotations.stream()
                    .anyMatch(a -> !a.HGVSp.isEmpty() && whitelist.get(a.Transcript).contains(a.HGVSp));

            final Predicate<VariantModel> predicate = v -> inPanel.test(v) ? !inBlacklist.test(v) : inWhitelist.test(v);
            final Predicate<GeneCopyNumber> copyNumberPredicate =
                    cnv -> cnvPredicates.stream().anyMatch(p -> p.test(cnv)) && cnv.minCopyNumber() < MAX_COPY_NUMBER_FOR_LOSS;
            final Predicate<HmfGenomeRegion> disruptionPredicate =
                    disruption -> disruptionPredicates.stream().anyMatch(p -> p.test(disruption));

            result.programs.put(program.getName(),
                    new BachelorProgram(program.getName(), predicate, copyNumberPredicate, disruptionPredicate));
        }

        if (iterateAllGenes) {
            result.querySet.putAll(allTranscriptsMap);
        }

        return result;
    }

    @NotNull
    private Collection<EligibilityReport> processVariant(final VariantContext variant, final String patient, final String sample,
            final EligibilityReport.ReportType type) {

        if (variant.isFiltered()) {
            return Collections.emptyList();
        }

        // we will skip when an ALT is not present in the sample
        final Genotype genotype = variant.getGenotype(sample);
        if (genotype == null || !(genotype.isHomVar() || genotype.isHet())) {
            return Collections.emptyList();
        }

        // TODO: do we need to verify specific ALTS have specific SnpEff effects

        final VariantModel model = VariantModel.from(variant);

        final List<String> matchingPrograms = programs.entrySet()
                .stream()
                .filter(program -> program.getValue().vcfProcessor.test(model))
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());

        // TODO: only add matching annotations of variant
        final String alts = variant.getAlternateAlleles().stream().map(Object::toString).collect(Collectors.joining("|"));
        final String effects = String.join("|", model.Annotations.stream().flatMap(a -> a.Effects.stream()).collect(Collectors.toSet()));
        final String genes = String.join("|", model.Annotations.stream().map(a -> a.GeneName).collect(Collectors.toSet()));

        return matchingPrograms.stream()
                .map(p -> ImmutableEligibilityReport.builder()
                        .patient(patient)
                        .source(type)
                        .program(p)
                        .id(variant.getID())
                        .genes(genes)
                        .chrom(variant.getContig())
                        .pos(variant.getStart())
                        .ref(variant.getReference().toString())
                        .alts(alts)
                        .effects(effects)
                        .build())
                .collect(Collectors.toList());
    }

    @NotNull
    Collection<EligibilityReport> processVCF(final String patient, final String sample, final EligibilityReport.ReportType type,
            final VCFFileReader reader) {

        final List<EligibilityReport> results = Lists.newArrayList();

        for (final HmfGenomeRegion region : querySet.values()) {
            final CloseableIterator<VariantContext> query =
                    reader.query(region.chromosome(), (int) region.geneStart(), (int) region.geneEnd());
            while (query.hasNext()) {
                final VariantContext variant = query.next();
                results.addAll(processVariant(variant, patient, sample, type));
            }
            query.close();
        }

        return results;
    }

    @NotNull
    public Collection<EligibilityReport> processCopyNumbers(final String patient, final List<GeneCopyNumber> copyNumbers) {

        final List<EligibilityReport> results = Lists.newArrayList();
        for (final GeneCopyNumber copyNumber : copyNumbers) {

            // TODO: verify the germline check
            final boolean isGermline = copyNumber.germlineHet2HomRegions() + copyNumber.germlineHomRegions() > 0;
            final List<String> matchingPrograms = programs.entrySet()
                    .stream()
                    .filter(program -> program.getValue().copyNumberProcessor.test(copyNumber))
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toList());

            final List<EligibilityReport> interimResults = matchingPrograms.stream()
                    .map(p -> ImmutableEligibilityReport.builder()
                            .patient(patient)
                            .source(isGermline ? GERMLINE_DELETION : SOMATIC_DELETION)
                            .program(p)
                            .id("")
                            .genes(copyNumber.gene())
                            .chrom(copyNumber.chromosome())
                            .pos(copyNumber.start())
                            .ref("")
                            .alts("")
                            .effects("")
                            .build())
                    .collect(Collectors.toList());

            results.addAll(interimResults);
        }

        return results;
    }

    private static int intron(final List<HmfExonRegion> exome, final GenomePosition position) {
        for (int i = 0; i < exome.size() - 1; i++) {
            if (position.position() > exome.get(i).end() && position.position() < exome.get(i + 1).start()) {
                return i;
            }
        }
        return -1;
    }

    private Collection<EligibilityReport> processStructuralVariant(final String patient, final GenomePosition position,
            final GenomePosition other, final StructuralVariantType svType) {

        final List<EligibilityReport> results = Lists.newArrayList();

        // TODO: can we do better than this performance wise? new map?
        for (final HmfGenomeRegion region : allGenesByChromosomeMap.get(position.chromosome())) {

            if (!region.contains(position)) {
                continue;
            }

            // skip non-inversion intronic variants
            if (region.contains(other) && svType != StructuralVariantType.INV) {
                final int intronStart = intron(region.exome(), position);
                final int intronEnd = intron(region.exome(), other);

                // the variant is intronic in a gene -- we will filter it
                if (intronStart >= 0 && intronStart == intronEnd) {
                    continue;
                }
            }

            programs.values()
                    .stream()
                    .filter(p -> p.disruptionProcessor.test(region))
                    .map(p -> ImmutableEligibilityReport.builder()
                            .patient(patient)
                            .source(SOMATIC_DISRUPTION)
                            .program(p.name)
                            .id("")
                            .genes(region.gene())
                            .chrom(region.chromosome())
                            .pos(position.position())
                            .ref("")
                            .alts("")
                            .effects("")
                            .build())
                    .forEach(results::add);

        }

        return results;
    }

    private Stream<EligibilityReport> processStructuralVariant(final String patient, final StructuralVariant structuralVariant) {
        final GenomePosition start = GenomePositions.create(structuralVariant.chromosome(true), structuralVariant.position(true));
        final GenomePosition end = GenomePositions.create(structuralVariant.chromosome(false), structuralVariant.position(false));

        final List<EligibilityReport> results = Lists.newArrayList();
        results.addAll(processStructuralVariant(patient, start, end, structuralVariant.type()));
        results.addAll(processStructuralVariant(patient, end, start, structuralVariant.type()));
        return results.stream();
    }

    @NotNull
    public Collection<EligibilityReport> processStructuralVariants(final String patient, final List<StructuralVariant> structuralVariants) {
        return structuralVariants.stream().flatMap(sv -> processStructuralVariant(patient, sv)).collect(Collectors.toList());
    }
}
