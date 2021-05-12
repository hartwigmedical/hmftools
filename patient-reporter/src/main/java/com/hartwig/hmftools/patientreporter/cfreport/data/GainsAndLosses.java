package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.protect.cnchromosome.CnPerChromosome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GainsAndLosses {

    private GainsAndLosses() {
    }

    private static final Logger LOGGER = LogManager.getLogger(GainsAndLosses.class);

    @NotNull
    public static List<ReportableGainLoss> sort(@NotNull List<ReportableGainLoss> reportableGainsAndLosses) {
        return reportableGainsAndLosses.stream().sorted((gainLoss1, gainLoss2) -> {
            String location1 = GeneUtil.zeroPrefixed(gainLoss1.chromosome() + gainLoss1.chromosomeBand());
            String location2 = GeneUtil.zeroPrefixed(gainLoss2.chromosome() + gainLoss2.chromosomeBand());

            if (location1.equals(location2)) {
                return gainLoss1.gene().compareTo(gainLoss2.gene());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static Set<String> amplifiedGenes(@NotNull List<ReportableGainLoss> reportableGainLosses) {
        Set<String> genes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainLosses) {
            if (gainLoss.interpretation() == CopyNumberInterpretation.FULL_GAIN
                    || gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_GAIN) {
                genes.add(gainLoss.gene());
            }
        }
        return genes;
    }

    @NotNull
    public static Set<String> lostGenes(@NotNull List<ReportableGainLoss> reportableGainLosses) {
        Set<String> genes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainLosses) {
            if (gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS
                    || gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS) {
                genes.add(gainLoss.gene());
            }
        }
        return genes;
    }

    public static double copyChromosomeArm(@Nullable CnPerChromosome cnPerChromosome, @NotNull String chromosome,
            @NotNull String chromosomeBand) {
        assert  cnPerChromosome != null;

        if (chromosome.equals("1") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr1p();
        } else if (chromosome.equals("1") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr1q();
        } else if (chromosome.equals("2") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr2p();
        } else if (chromosome.equals("2") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr2q();
        } else if (chromosome.equals("3") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr3p();
        } else if (chromosome.equals("3") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr3q();
        } else if (chromosome.equals("4") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr4p();
        } else if (chromosome.equals("4") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr4q();
        } else if (chromosome.equals("5") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr5p();
        } else if (chromosome.equals("5") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr5q();
        } else if (chromosome.equals("6") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr6p();
        } else if (chromosome.equals("6") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr6q();
        } else if (chromosome.equals("7") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr7p();
        } else if (chromosome.equals("7") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr7q();
        } else if (chromosome.equals("8") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr8p();
        } else if (chromosome.equals("8") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr8q();
        } else if (chromosome.equals("9") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr9p();
        } else if (chromosome.equals("9") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr9q();
        } else if (chromosome.equals("10") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr10p();
        } else if (chromosome.equals("10") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr10q();
        } else if (chromosome.equals("11") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr11p();
        } else if (chromosome.equals("11") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr11q();
        } else if (chromosome.equals("12") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr12p();
        } else if (chromosome.equals("12") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr12q();
        } else if (chromosome.equals("13") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr13p();
        } else if (chromosome.equals("13") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr13q();
        } else if (chromosome.equals("14") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr14p();
        } else if (chromosome.equals("14") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr14q();
        } else if (chromosome.equals("15") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr15p();
        } else if (chromosome.equals("15") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr15q();
        } else if (chromosome.equals("16") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr16p();
        } else if (chromosome.equals("16") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr16q();
        } else if (chromosome.equals("17") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr17p();
        } else if (chromosome.equals("17") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr17q();
        } else if (chromosome.equals("18") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr18p();
        } else if (chromosome.equals("18") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr18q();
        } else if (chromosome.equals("19") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr19p();
        } else if (chromosome.equals("19") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr19q();
        } else if (chromosome.equals("20") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr20p();
        } else if (chromosome.equals("20") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr20q();
        } else if (chromosome.equals("21") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr21p();
        } else if (chromosome.equals("21") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr21q();
        } else if (chromosome.equals("22") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chr22p();
        } else if (chromosome.equals("22") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chr22q();
        } else if (chromosome.equals("X") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chrXp();
        } else if (chromosome.equals("X") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chrXq();
        } else if (chromosome.equals("Y") && chromosomeBand.contains("p")) {
            return cnPerChromosome.chrYp();
        } else if (chromosome.equals("Y") && chromosomeBand.contains("q")) {
            return cnPerChromosome.chrYq();
        } else {
            LOGGER.warn("Chromosome copy number could not be extracted");
            return 0;
        }

    }
}
