package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.IOException;
import java.util.List;

import org.apache.commons.compress.utils.Lists;

public class DriverGenePanelConversion {

    public static void main(String[] args) throws IOException {
        String inputFile = "/Users/jon/hmf/resources/DriverGenePanel.hg19.tsv";
        String outputFile = "/Users/jon/hmf/resources/DriverGenePanel.hg38.tsv";

        DndsGeneNameMap geneNameMap = new DndsGeneNameMap();
        List<DriverGene> inputDriverGenes = DriverGeneFile.read(inputFile);
        List<DriverGene> outputDriverGenes = Lists.newArrayList();
        for (DriverGene input : inputDriverGenes) {
            final String hg19Gene = input.gene();
            final String hg38Gene = geneNameMap.hg38Gene(input.gene());
            if (input.reportVariant() && !geneNameMap.isValidHg19Gene(hg19Gene)) {
                System.out.println("Bad gene: " + hg19Gene);
            }

            if (hg19Gene.equals("LINC00290") || hg19Gene.equals("LINC01001")) {
                outputDriverGenes.add(input);
            } else if (hg38Gene.equals("NA")) {
                System.out.println("Excluding: " + hg19Gene);
            } else {
                DriverGene converted = ImmutableDriverGene.builder().from(input).gene(hg38Gene).build();
                outputDriverGenes.add(converted);
            }
        }
        DriverGeneFile.write(outputFile, outputDriverGenes);
    }
}
