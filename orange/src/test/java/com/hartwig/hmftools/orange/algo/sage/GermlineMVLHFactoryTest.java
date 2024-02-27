package com.hartwig.hmftools.orange.algo.sage;

import static org.apache.commons.math3.util.Precision.EPSILON;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;

import org.junit.Test;

public class GermlineMVLHFactoryTest
{
    @Test
    public void parseGermlineMVLHCorrectly()
    {
        String reliableMVLHGene = "FAKE1";
        String unReliableMVLHGene = "FAKE2";

        List<String> lines = List.of(
                "gene\tchromosome\tposStart\tposEnd\tmissedVariantLikelihood\tDR_0\tDR_1\tDR_2\tDR_3\tDR_4\tDR_5\tDR_6\tDR_7\tDR_8\tDR_9\tDR_10\tDR_11\tDR_12\tDR_13\tDR_14\tDR_15\tDR_16\tDR_17\tDR_18\tDR_19\tDR_20\tDR_21\tDR_22\tDR_23\tDR_24\tDR_25\tDR_26\tDR_27\tDR_28\tDR_29\tDR_30_39\tDR_40_49\tDR_50_59\tDR_60_69\tDR_70_79\tDR_80_89\tDR_90_99\tDR_100_149\tDR_150_199\tDR_200_249\tDR_250_299\tDR_300_349\tDR_350_399\tDR_400_449\tDR_450_499\tDR_500_599\tDR_600_699\tDR_700_799\tDR_800_899\tDR_900_999\tDR_1000_1099\tDR_1100_1199\tDR_1200_1299\tDR_1300_1399\tDR_1400_1499\tDR_1500_1599\tDR_1600_1699\tDR_1700_1799\tDR_1800_1899\tDR_1900_1999\tDR_2000_2999\tDR_3000_3999\tDR_4000_4999\tDR_5000_5999\tDR_6000_6999\tDR_7000_7999\tDR_8000_8999\tDR_9000_9999\tDR_10000",
                reliableMVLHGene + "\t7\t87133559\t87229500\t0.0002\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t9\t41\t35\t48\t41\t85\t76\t132\t132\t162\t222\t237\t2689\t454\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0",
                unReliableMVLHGene +"\t1\t120436587\t120438959\t1.0000\t2373\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"
        );

        DriverGene driverGene1 = ImmutableDriverGene.builder()
                .gene(reliableMVLHGene)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportDisruption(true)
                .reportAmplification(false)
                .reportSomaticHotspot(true)
                .likelihoodType(DriverCategory.TSG)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .reportGermlineDeletion(DriverGeneGermlineReporting.NONE)
                .reportGermlineDisruption(DriverGeneGermlineReporting.NONE)
                .reportPGX(false)
                .build();
        DriverGene driverGene2 = ImmutableDriverGene.builder()
                .gene(unReliableMVLHGene)
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(false)
                .reportAmplification(true)
                .reportSomaticHotspot(false)
                .likelihoodType(DriverCategory.ONCO)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .reportGermlineDeletion(DriverGeneGermlineReporting.NONE)
                .reportGermlineDisruption(DriverGeneGermlineReporting.NONE)
                .reportPGX(false)
                .build();
        List<DriverGene> driverGenes = List.of(driverGene1, driverGene2);

        Map<String, Double> mvlhPerGene = GermlineMVLHFactory.parseMVLHPerGene(lines, driverGenes);

        assertNotNull(mvlhPerGene);
        assertEquals(0.0002D, mvlhPerGene.get(reliableMVLHGene), EPSILON);
        assertFalse(mvlhPerGene.containsKey(unReliableMVLHGene));
    }
}