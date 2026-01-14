package com.hartwig.hmftools.datamodel.finding;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.purple.PurpleCopyNumber;

import org.jetbrains.annotations.NotNull;

final class ChromosomeArmCopyNumberMap
{
    private enum ChromosomeArm
    {
        P_ARM,
        Q_ARM,
    }

    private record CopyNumberKey(String chromosome, ChromosomeArm arm) {
    }

    private final Map<CopyNumberKey, Double> cnPerChromosomeArm;

    public static ChromosomeArmCopyNumberMap create(@NotNull Iterable<PurpleCopyNumber> copyNumbers,
            @NotNull final OrangeRefGenomeVersion refGenomeVersion) {
        return new ChromosomeArmCopyNumberMap(extractCnPerChromosomeArm(copyNumbers, refGenomeVersion));
    }

    private ChromosomeArmCopyNumberMap(Map<CopyNumberKey, Double> cnPerChromosomeArm)
    {
        this.cnPerChromosomeArm = cnPerChromosomeArm;
    }

    public double chromosomeArmCopyNumber(@NotNull String chromosome, @NotNull String chromosomeBand) {
        Double copyNumber = cnPerChromosomeArm.get(new CopyNumberKey(chromosome, getChromosomeArm(chromosomeBand)));
        if (copyNumber == null) {
            throw new IllegalArgumentException("Copy number not found for chromosome band: " + chromosomeBand + "!");
        }
        return copyNumber;
    }

    @NotNull
    static ChromosomeArm getChromosomeArm(@NotNull String chromosomeBand) {
        ChromosomeArm chromosomeArm;
        if (chromosomeBand.startsWith("p")) {
            chromosomeArm = ChromosomeArm.P_ARM;
        } else if (chromosomeBand.startsWith("q")) {
            chromosomeArm = ChromosomeArm.Q_ARM;
        } else {
            throw new IllegalArgumentException("Chromosome arm could not be resolved from band: " + chromosomeBand + "!");
        }
        return chromosomeArm;
    }

    @NotNull
    private static Map<CopyNumberKey, Double> extractCnPerChromosomeArm(@NotNull Iterable<PurpleCopyNumber> copyNumbers,
            @NotNull final OrangeRefGenomeVersion refGenomeVersion)
    {
        RefGenomeCoordinates refGenomeCoordinates = RefGenomeCoordinates.refGenomeCoordinates(refGenomeVersion);

        Map<CopyNumberKey, Double> cnPerChromosomeArmData = new HashMap<>();
        for(String chromosome : refGenomeCoordinates.Lengths.keySet())
        {
            Map<ChromosomeArm, GenomeRegion> genomeRegion = determineArmRegions(chromosome, refGenomeCoordinates);

            for(ChromosomeArm arm : ChromosomeArm.values())
            {
                double copyNumberArm = 0;

                GenomeRegion armGenomeRegion = genomeRegion.get(arm);

                if(armGenomeRegion != null)
                {
                    for (PurpleCopyNumber purpleCopyNumber : copyNumbers) {
                        String copyNumberChromosome = purpleCopyNumber.chromosome();

                        if (copyNumberChromosome.equals(chromosome) && overlaps(armGenomeRegion, purpleCopyNumber)) {
                            double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                            int totalLengthSegment = bases(purpleCopyNumber);
                            copyNumberArm += (copyNumber * totalLengthSegment) / bases(armGenomeRegion);
                        }
                    }
                }

                copyNumberArm = Doubles.round(copyNumberArm, 2);

                cnPerChromosomeArmData.put(new CopyNumberKey(chromosome, arm), copyNumberArm);
            }
        }
        return cnPerChromosomeArmData;
    }

    @NotNull
    private static Map<ChromosomeArm, GenomeRegion> determineArmRegions(@NotNull String chromosome,
            @NotNull RefGenomeCoordinates refGenomeCoordinates)
    {
        int centromerePos = refGenomeCoordinates.centromeres().get(chromosome);
        int chrLength = refGenomeCoordinates.lengths().get(chromosome);

        Map<ChromosomeArm, GenomeRegion> chromosomeArmGenomeRegionMap = new HashMap<>();

        GenomeRegion partBeforeCentromere = new GenomeRegion(chromosome, 1, centromerePos);
        GenomeRegion partAfterCentromere = new GenomeRegion(chromosome, centromerePos + 1, chrLength);

        if (bases(partBeforeCentromere) < bases(partAfterCentromere))
        {
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.P_ARM, partBeforeCentromere);
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.Q_ARM, partAfterCentromere);
        } else
        {
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.P_ARM, partAfterCentromere);
            chromosomeArmGenomeRegionMap.put(ChromosomeArm.Q_ARM, partBeforeCentromere);
        }
        return chromosomeArmGenomeRegionMap;
    }

    private static boolean overlaps(@NotNull GenomeRegion region, @NotNull PurpleCopyNumber purpleCopyNumber)
    {
        return purpleCopyNumber.chromosome().equals(region.chromosome()) && purpleCopyNumber.end() > region.start()
                && purpleCopyNumber.start() < region.end();
    }

    private static int bases(@NotNull GenomeRegion region)
    {
        return bases(region.start(), region.end());
    }

    private static int bases(@NotNull PurpleCopyNumber purpleCopyNumber)
    {
        return bases(purpleCopyNumber.start(), purpleCopyNumber.end());
    }

    private static int bases(int start, int end)
    {
        return 1 + end - start;
    }

    public enum RefGenomeCoordinates
    {
        COORDS_37(v37Lengths(), v37centromeres()),
        COORDS_38(v38Lengths(), v38centromeres());

        public final Map<String, Integer> Lengths;
        public final Map<String, Integer> Centromeres;

        public static RefGenomeCoordinates refGenomeCoordinates(final OrangeRefGenomeVersion refGenomeVersion)
        {
            return switch (refGenomeVersion) {
                case V37 -> COORDS_37;
                case V38 -> COORDS_38;
            };
        }

        RefGenomeCoordinates(@NotNull final Map<String, Integer> lengths, @NotNull final Map<String ,Integer> centromeres)
        {
            Lengths = lengths;
            Centromeres = centromeres;
        }

        @NotNull
        public Map<String, Integer> lengths() {
            return Lengths;
        }

        @NotNull
        public Map<String, Integer> centromeres() {
            return Centromeres;
        }

        private static Map<String, Integer> v37Lengths()
        {
            Map<String, Integer> lengths = new HashMap<>();
            lengths.put("1", 249250621);
            lengths.put("2", 243199373);
            lengths.put("3", 198022430);
            lengths.put("4", 191154276);
            lengths.put("5", 180915260);
            lengths.put("6", 171115067);
            lengths.put("7", 159138663);
            lengths.put("8", 146364022);
            lengths.put("9", 141213431);
            lengths.put("10", 135534747);
            lengths.put("11", 135006516);
            lengths.put("12", 133851895);
            lengths.put("13", 115169878);
            lengths.put("14", 107349540);
            lengths.put("15", 102531392);
            lengths.put("16", 90354753);
            lengths.put("17", 81195210);
            lengths.put("18", 78077248);
            lengths.put("19", 59128983);
            lengths.put("20", 63025520);
            lengths.put("21", 48129895);
            lengths.put("22", 51304566);
            lengths.put("X", 155270560);
            lengths.put("Y", 59373566);
            return lengths;
        }

        private static Map<String, Integer> v38Lengths()
        {
            Map<String, Integer> lengths = new HashMap<>();
            lengths.put("chr1", 248956422);
            lengths.put("chr2", 242193529);
            lengths.put("chr3", 198295559);
            lengths.put("chr4", 190214555);
            lengths.put("chr5", 181538259);
            lengths.put("chr6", 170805979);
            lengths.put("chr7", 159345973);
            lengths.put("chr8", 145138636);
            lengths.put("chr9", 138394717);
            lengths.put("chr10", 133797422);
            lengths.put("chr11", 135086622);
            lengths.put("chr12", 133275309);
            lengths.put("chr13", 114364328);
            lengths.put("chr14", 107043718);
            lengths.put("chr15", 101991189);
            lengths.put("chr16", 90338345);
            lengths.put("chr17", 83257441);
            lengths.put("chr18", 80373285);
            lengths.put("chr19", 58617616);
            lengths.put("chr20", 64444167);
            lengths.put("chr21", 46709983);
            lengths.put("chr22", 50818468);
            lengths.put("chrX", 156040895);
            lengths.put("chrY", 57227415);
            return lengths;
        }

        private static Map<String, Integer> v37centromeres()
        {
            Map<String, Integer> lengths = new HashMap<>();
            lengths.put("1", 123035434);
            lengths.put("2", 93826171);
            lengths.put("3", 92004854);
            lengths.put("4", 51160117);
            lengths.put("5", 47905641);
            lengths.put("6", 60330166);
            lengths.put("7", 59554331);
            lengths.put("8", 45338887);
            lengths.put("9", 48867679);
            lengths.put("10", 40754935);
            lengths.put("11", 53144205);
            lengths.put("12", 36356694);
            lengths.put("13", 17500000);
            lengths.put("14", 17500000);
            lengths.put("15", 18500000);
            lengths.put("16", 36835801);
            lengths.put("17", 23763006);
            lengths.put("18", 16960898);
            lengths.put("19", 26181782);
            lengths.put("20", 27869569);
            lengths.put("21", 12788129);
            lengths.put("22", 14500000);
            lengths.put("X", 60132012);
            lengths.put("Y", 11604553);
            return lengths;
        }

        private static Map<String, Integer> v38centromeres()
        {
            Map<String, Integer> lengths = new HashMap<>();
            lengths.put("chr1", 123605523);
            lengths.put("chr2", 93139351);
            lengths.put("chr3", 92214016);
            lengths.put("chr4", 50726026);
            lengths.put("chr5", 48272854);
            lengths.put("chr6", 59191911);
            lengths.put("chr7", 59498944);
            lengths.put("chr8", 44955505);
            lengths.put("chr9", 44377363);
            lengths.put("chr10", 40640102);
            lengths.put("chr11", 52751711);
            lengths.put("chr12", 35977330);
            lengths.put("chr13", 17025624);
            lengths.put("chr14", 17086762);
            lengths.put("chr15", 18362627);
            lengths.put("chr16", 37295920);
            lengths.put("chr17", 24849830);
            lengths.put("chr18", 18161053);
            lengths.put("chr19", 25844927);
            lengths.put("chr20", 28237290);
            lengths.put("chr21", 11890184);
            lengths.put("chr22", 14004553);
            lengths.put("chrX", 60509061);
            lengths.put("chrY", 10430492);
            return lengths;
        }
    }

    private record GenomeRegion(String chromosome, int start, int end) {}
}
