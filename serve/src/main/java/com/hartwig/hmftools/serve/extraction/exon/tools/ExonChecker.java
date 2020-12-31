package com.hartwig.hmftools.serve.extraction.exon.tools;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExonFile;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

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

        List<HmfTranscriptRegion> transcripts = HmfGenePanelSupplier.allGeneList37();
        LOGGER.info("The size of the file is {}", exons.size());

        String chromosome = null;
        Long start = null;
        Long end = null;
        String exonEnsemblId = null;
        String gene = null;
        for (KnownExon exon : exons) {
            chromosome = exon.annotation().chromosome();
            start = exon.annotation().start();
            end = exon.annotation().end();
            exonEnsemblId = exon.annotation().exonEnsemblId();
            gene = exon.annotation().gene();
            List<HmfTranscriptRegion> transriptsGenes = Lists.newArrayList();
            List<String> transriptsExonIds = Lists.newArrayList();

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
                        LOGGER.debug("The exon Id of SERVE {}, is not the same as in GRch 37{} but is placed in other loop",
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
        LOGGER.info("Checking exons");

        LOGGER.info("Done!");

    }

}
