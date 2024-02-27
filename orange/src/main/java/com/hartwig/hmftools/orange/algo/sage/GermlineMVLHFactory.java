package com.hartwig.hmftools.orange.algo.sage;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.sage.GeneDepthFile;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GermlineMVLHFactory
{
    @NotNull
    public static Map<String, Double> loadGermlineMVLHPerGene(@NotNull String sageGermlineGeneCoverageTsv,
            @NotNull List<DriverGene> driverGenes) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(sageGermlineGeneCoverageTsv).toPath());
        return parseMVLHPerGene(lines, driverGenes);
    }

    @NotNull
    @VisibleForTesting
    static Map<String, Double> parseMVLHPerGene(@NotNull List<String> lines, @NotNull List<DriverGene> driverGenes)
    {
        String header = lines.get(0);

        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        int geneIndex = fieldsIndexMap.get(GeneDepthFile.COL_GENE);
        int mvlhIndex = fieldsIndexMap.get(GeneDepthFile.COL_MV_LIKELIHOOD);

        Map<String, Double> mvlhPerGene = Maps.newTreeMap();
        for(String line : lines.subList(1, lines.size()))
        {
            String[] values = line.split(TSV_DELIM);
            String gene = values[geneIndex];
            double mvlh = Double.parseDouble(values[mvlhIndex]);

            DriverGene matchingDriverGene = findMatchingDriverGene(gene, driverGenes);
            if(matchingDriverGene != null && hasReliableMVLH(matchingDriverGene))
            {
                mvlhPerGene.put(gene, mvlh);
            }
        }
        return mvlhPerGene;
    }

    @Nullable
    private static DriverGene findMatchingDriverGene(@NotNull String geneName, @NotNull List<DriverGene> driverGenes)
    {
        List<DriverGene> matchingDriverGenes = driverGenes.stream().filter(d -> d.gene().equals(geneName)).collect(Collectors.toList());
        if(matchingDriverGenes.size() == 1)
        {
            return matchingDriverGenes.get(0);
        }
        else if(matchingDriverGenes.isEmpty())
        {
            return null;
        }
        else
        {
            throw new IllegalStateException("More than one driver gene found with name " + geneName);
        }
    }

    private static boolean hasReliableMVLH(@NotNull DriverGene driverGene)
    {
        // Note: fix this when SAGE germline starts running genome wide
        return driverGene.reportSomatic() || driverGene.reportGermline() || driverGene.reportPGX();
    }
}
