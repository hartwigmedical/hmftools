package com.hartwig.hmftools.serve.extraction.exon.tools;

import java.io.IOException;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExonFile;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ExonChecker {
    private static final Logger LOGGER = LogManager.getLogger(ExonChecker.class);
    private static final boolean LOG_DEBUG = true;

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE exon checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String knownExonsTsv = System.getProperty("user.home") + "/hmf/tmp/serve/KnownExons.SERVE.37.tsv";
        List<KnownExon> exons = KnownExonFile.read(knownExonsTsv);
        Random random = new Random();
        List<HmfTranscriptRegion> transcripts = HmfGenePanelSupplier.allGeneList37();
        LOGGER.info("The size of the file is {}", exons.size());

        List<String> bases = Lists.newArrayList();
        String randomAltBase = null;
        for (KnownExon exon : exons) {
            String chromosome = exon.annotation().chromosome();
            Long start = exon.annotation().start();
            Long end = exon.annotation().end();
            String exonEnsemblId = exon.annotation().exonEnsemblId();
            String gene = exon.annotation().gene();
            List<HmfTranscriptRegion> transriptsGenes = Lists.newArrayList();
            List<String> transriptsExonIds = Lists.newArrayList();
            MutationTypeFilter mutationTypeFilter = exon.annotation().mutationType();

            //Extract genomicpositions with SNPEFF
            Long l = new Long(end - start + 1);
            Long bewtweenNumber = start + random.nextInt(l.intValue());
            LOGGER.info(bewtweenNumber);
            List<Long> genomicPositions = Lists.newArrayList(start, bewtweenNumber, end);

            for (Long genomicPosition : genomicPositions) {
                String extractRefBaseOfPosition = extractRefBaseOfGenomicPosition(chromosome, genomicPosition);

                if (extractRefBaseOfPosition.equals("A")) {
                    bases = Lists.newArrayList("C", "T", "G");
                    randomAltBase = bases.get(random.nextInt(3));
                } else if (extractRefBaseOfPosition.equals("C")) {
                    bases = Lists.newArrayList("A", "T", "G");
                    randomAltBase = bases.get(random.nextInt(3));
                } else if (extractRefBaseOfPosition.equals("T")) {
                    bases = Lists.newArrayList("C", "A", "G");
                    randomAltBase = bases.get(random.nextInt(3));
                } else if (extractRefBaseOfPosition.equals("G")) {
                    bases = Lists.newArrayList("C", "T", "A");
                    randomAltBase = bases.get(random.nextInt(3));
                }

                extactMutationType(extractRefBaseOfPosition, randomAltBase, chromosome, genomicPosition, mutationTypeFilter);
            }

            //Extract genomic positions of GRCH 37
            for (HmfTranscriptRegion region : transcripts) {
                if (region.gene().equals(gene)) {
                    transriptsGenes.add(region);
                    for (HmfExonRegion exonRegion : region.exome()) {
                        transriptsExonIds.add(exonRegion.exonID());
                    }
                }
            }

            for (HmfTranscriptRegion transriptsGene : transriptsGenes) {
                for (HmfExonRegion exonRegion : transriptsGene.exome()) {
                    if (transriptsExonIds.contains(exonRegion.exonID())) {
                        LOGGER.debug("The exon Id of SERVE {}, is not the same as in GRch 37{} but exon ID will be used",
                                exonEnsemblId,
                                exonRegion.exonID());
                    } else if (exonRegion.exonID().equals(exonEnsemblId)) {
                        Long exonStart = exonRegion.start() - 5;
                        Long exonEnd = exonRegion.end() + 5;
                        String exonChromosome = exonRegion.chromosome();
                        if (!exonStart.equals(start)) {
                            LOGGER.warn("The exon start postion of SERVE {} is not the same as in GRch 37 {} on exon ID {}",
                                    start,
                                    exonStart,
                                    exonEnsemblId);
                        }
                        if (!exonEnd.equals(end)) {
                            LOGGER.warn("The exon end postion of SERVE {} is not the same as in GRch 37 {} on exon ID {}",
                                    end,
                                    exonEnd,
                                    exonEnsemblId);
                        }
                        if (!exonChromosome.equals(chromosome)) {
                            LOGGER.warn("The exon chromosome of SERVE {} is not the same as in GRch 37 {} on exon ID {}",
                                    chromosome,
                                    exonChromosome,
                                    exonEnsemblId);
                        }
                    } else {
                        LOGGER.warn("The exon Id of SERVE {}, is not the same as in GRch 37{}", exonEnsemblId, exonRegion.exonID());
                    }
                }
            }
        }
        LOGGER.info("All exons are checked!");

        LOGGER.info("Done!");
    }

    private static String extractRefBaseOfGenomicPosition(@Nullable String chromosome, Long genomicPosition) {
        String genomicPositionBase = chromosome + ":" + genomicPosition + "-" + genomicPosition;
        //TODO extract base with genomic posiiton
        return genomicPositionBase;
    }

    private static void extactMutationType(@Nullable String extractRefBaseOfPosition, @Nullable String randomAltBase,
            @Nullable String chromosome, Long position, @NotNull MutationTypeFilter mutationTypeFilter) {
        //TODO extract codon annotation
        //TODO extract mutation type filter
    }
}
