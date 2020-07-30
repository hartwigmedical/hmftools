package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihoodSupplier;

import org.jetbrains.annotations.NotNull;

public class DriverGeneFactory {

    @NotNull
    private static DriverGene create(@NotNull final String gene) {
        return ImmutableDriverGene.builder()
                .gene(gene)
                .deletionBand(null)
                .reportMissense(true)
                .reportTruncation(false)
                .reportSplice(false)
                .reportDisruption(false)
                .reportAmplification(false)
                .favorMultiHitAndBiallelic(false)
                .build();
    }

    @NotNull
    private static Set<String> amplificationTargets() {
        final InputStream inputStream = DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream("/cna/AmplificationTargets.tsv");
        return new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toSet());
    }

    @NotNull
    private static Map<String, String> deletionTargets() {
        final Map<String, String> result = Maps.newHashMap();
        final InputStream inputStream = DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream("/cna/DeletionTargets.tsv");
        new BufferedReader(new InputStreamReader(inputStream)).lines().forEach(line -> {
            final String[] values = line.split("\t");
            result.put(values[0], values[1]);
        });

        return result;
    }

    public static void main(String[] args) throws IOException {

        final Map<String, DriverGene> all = Maps.newHashMap();
        final Set<String> oncoGenes = DndsDriverGeneLikelihoodSupplier.oncoLikelihood().keySet();
        for (String gene : oncoGenes) {
            DriverGene driver = all.computeIfAbsent(gene, DriverGeneFactory::create);
            DriverGene updated = ImmutableDriverGene.builder().from(driver).reportMissense(true).reportAmplification(true).build();
            all.put(gene, updated);
        }

        final Set<String> tsGenes = DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet();
        for (String gene : tsGenes) {
            DriverGene driver = all.computeIfAbsent(gene, DriverGeneFactory::create);
            DriverGene updated = ImmutableDriverGene.builder()
                    .from(driver)
                    .reportMissense(true)
                    .reportTruncation(true)
                    .reportSplice(true)
                    .reportDisruption(true)
                    .favorMultiHitAndBiallelic(true)
                    .build();
            all.put(gene, updated);
        }

        for (String gene : amplificationTargets()) {
            DriverGene driver = all.computeIfAbsent(gene, DriverGeneFactory::create);
            DriverGene updated = ImmutableDriverGene.builder().from(driver).reportAmplification(true).build();
            all.put(gene, updated);
        }

        Map<String, String> deletionTargets = deletionTargets();
        for (Map.Entry<String, String> entry : deletionTargets.entrySet()) {
            DriverGene driver = all.computeIfAbsent(entry.getKey(), DriverGeneFactory::create);
            DriverGene updated = ImmutableDriverGene.builder()
                    .from(driver)
                    .reportDisruption(true)
                    .deletionBand(entry.getValue().equals("NA") ? null : entry.getValue())
                    .build();
            all.put(entry.getKey(), updated);
        }

        List<DriverGene> resultList = new ArrayList<>(all.values());
        Collections.sort(resultList);

        DriverGeneFile.write("/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/drivercatalog/DriverGenePanel.tsv",
                new ArrayList<>(resultList));
    }

}
