package com.hartwig.hmftools.common.centromeres;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.refgenome.RefGenome;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.jetbrains.annotations.NotNull;

public class Centromeres {

    private final Map<Chromosome, Long> centromeres;

    public Centromeres(@NotNull final RefGenome refGenome) {
        switch (refGenome) {
            case HG38:
                centromeres = hg38();
                break;
            default:
                centromeres = hg19();
                break;

        }
    }

    public Centromeres(@NotNull final Map<Chromosome, Long> centromeres) {
        this.centromeres = centromeres;
    }

    public Map<Chromosome, Long> centromeres() {
        return centromeres;
    }

    @Deprecated
    @NotNull
    public static Map<String, GenomeRegion> grch37() {
        final Map<String, GenomeRegion> result = Maps.newHashMap();
        addCentromere(result, "1", 121535434, 124535434);
        addCentromere(result, "2", 92326171, 95326171);
        addCentromere(result, "3", 90504854, 93504854);
        addCentromere(result, "4", 49660117, 52660117);
        addCentromere(result, "5", 46405641, 49405641);
        addCentromere(result, "6", 58830166, 61830166);
        addCentromere(result, "7", 58054331, 61054331);
        addCentromere(result, "8", 43838887, 46838887);
        addCentromere(result, "9", 47367679, 50367679);
        addCentromere(result, "10", 39254935, 42254935);
        addCentromere(result, "11", 51644205, 54644205);
        addCentromere(result, "12", 34856694, 37856694);
        addCentromere(result, "13", 16000000, 19000000);
        addCentromere(result, "14", 16000000, 19000000);
        addCentromere(result, "15", 17000000, 20000000);
        addCentromere(result, "16", 35335801, 38335801);
        addCentromere(result, "17", 22263006, 25263006);
        addCentromere(result, "18", 15460898, 18460898);
        addCentromere(result, "19", 24681782, 27681782);
        addCentromere(result, "20", 26369569, 29369569);
        addCentromere(result, "21", 11288129, 14288129);
        addCentromere(result, "22", 13000000, 16000000);
        addCentromere(result, "X", 58632012, 61632012);
        addCentromere(result, "Y", 10104553, 13104553);

        return result;
    }

    private static void addCentromere(@NotNull final Map<String, GenomeRegion> map, @NotNull final String chromosome, final long start,
            final long end) {
        map.put(chromosome, GenomeRegionFactory.create(chromosome, start, end));
    }

    @NotNull
    public static Map<Chromosome, Long> hg19() {
        final Map<Chromosome, Long> result = Maps.newHashMap();

        result.put(HumanChromosome._1, 123035434L);
        result.put(HumanChromosome._2, 93826171L);
        result.put(HumanChromosome._3, 92004854L);
        result.put(HumanChromosome._4, 51160117L);
        result.put(HumanChromosome._5, 47905641L);
        result.put(HumanChromosome._6, 60330166L);
        result.put(HumanChromosome._7, 59554331L);
        result.put(HumanChromosome._8, 45338887L);
        result.put(HumanChromosome._9, 48867679L);
        result.put(HumanChromosome._10, 40754935L);
        result.put(HumanChromosome._11, 53144205L);
        result.put(HumanChromosome._12, 36356694L);
        result.put(HumanChromosome._13, 17500000L);
        result.put(HumanChromosome._14, 17500000L);
        result.put(HumanChromosome._15, 18500000L);
        result.put(HumanChromosome._16, 36835801L);
        result.put(HumanChromosome._17, 23763006L);
        result.put(HumanChromosome._18, 16960898L);
        result.put(HumanChromosome._19, 26181782L);
        result.put(HumanChromosome._20, 27869569L);
        result.put(HumanChromosome._21, 12788129L);
        result.put(HumanChromosome._22, 14500000L);
        result.put(HumanChromosome._X, 60132012L);
        result.put(HumanChromosome._Y, 11604553L);

        return result;
    }

    @NotNull
    private static Map<Chromosome, Long> hg38() {
        final Map<Chromosome, Long> result = Maps.newHashMap();

        result.put(HumanChromosome._1, 123605523L);
        result.put(HumanChromosome._2, 93139351L);
        result.put(HumanChromosome._3, 92214016L);
        result.put(HumanChromosome._4, 50726026L);
        result.put(HumanChromosome._5, 48272854L);
        result.put(HumanChromosome._6, 59191911L);
        result.put(HumanChromosome._7, 59498944L);
        result.put(HumanChromosome._8, 44955505L);
        result.put(HumanChromosome._9, 44377363L);
        result.put(HumanChromosome._10, 40640102L);
        result.put(HumanChromosome._11, 52751711L);
        result.put(HumanChromosome._12, 35977330L);
        result.put(HumanChromosome._13, 17025624L);
        result.put(HumanChromosome._14, 17086762L);
        result.put(HumanChromosome._15, 18362627L);
        result.put(HumanChromosome._16, 37295920L);
        result.put(HumanChromosome._17, 24849830L);
        result.put(HumanChromosome._18, 18161053L);
        result.put(HumanChromosome._19, 25844927L);
        result.put(HumanChromosome._20, 28237290L);
        result.put(HumanChromosome._21, 11890184L);
        result.put(HumanChromosome._22, 14004553L);
        result.put(HumanChromosome._X, 60509061L);
        result.put(HumanChromosome._Y, 10430492L);

        return result;
    }

}
