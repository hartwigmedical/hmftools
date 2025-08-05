package com.hartwig.hmftools.pavereverse;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class TestDataCreator
{
//    @Test
    public void createReducedEnsemblDataSet() throws IOException
    {
        File fullEnsemblDataDir = new File("/Users/timlavers/work/data/v6_0/ref/38/common/ensembl_data");
        File outputDir = new File("/Users/timlavers/work/code/hmftools/purple/src/test/resources/ensembl");

        final Set<String> geneNames = Set.of(
                "ETV6,"
//                "TET2,",
//                "VHL,",
//                "BRAF,",
//                "ADCK2,",
//                "EGFR,",
//                "RNU1-82P,",
//                "ZYX,",
//                "PIK3R1",
//                "ARID1A",
//                "KIT",
//                "BRCA1",
//                "DYRK1A",
//                "TATDN2",
//                "TERT",
//                "SRGAP2,",
//                "AP5Z1,"
        );
//        copyLinesMatching(new File(fullEnsemblDataDir, "ensembl_protein_features.csv"), outputDir, Set.of());
                Set<String> geneIds = Set.of(
                        "ENSG00000198793",
                        "ENSG00000117713",
                        "ENSG00000266028",
                        "ENSG00000242802",
                        "ENSG00000146648",
                        "ENSG00000133597",
                        "ENSG00000157764",
                        "ENSG00000212153",
                        "ENSG00000159840",
                        "ENSG00000071564",
        "ENSG00000139083"
        //                "ENSG00000168769", // TET2
        //                "ENSG00000134086",  // VHL
        //                "ENSG00000133597", // ADCK2
        //                "ENSG00000157764", // BRAF
        //                "ENSG00000146648", // EGFR
        //                "ENSG00000212153", // RNU-82P
        //                "ENSG00000159840", // ZYX
        //                "ENSG00000145675", // PIK3R1
        //                "ENSG00000117713", // ARID1A
        //                "ENSG00000157404", // KIT
        //                "ENSG00000012048",  // BRCA1
        //                "ENSG00000157540", // DYRK1A
        //                "ENSG00000157014", // TATDN2
        //                "ENSG00000164362", // TERT
        //                "ENSG00000163486", // SRGAP2
        //                "ENSG00000242802" // AP5Z1
                );
//        Set<String> geneIds = Set.of(
//                "ENSG00000198793", // MTOR
//                "ENSG00000168769", // TET2
//                "ENSG00000134086",  // VHL
//                "ENSG00000133597", // ADCK2
//                "ENSG00000157764", // BRAF
//                "ENSG00000146648", // EGFR
//                "ENSG00000212153", // RNU-82P
//                "ENSG00000159840", // ZYX
//                "ENSG00000145675", // PIK3R1
//                "ENSG00000117713", // ARID1A
//                "ENSG00000157404", // KIT
//                "ENSG00000012048",  // BRCA1
//                "ENSG00000157540", // DYRK1A
//                "ENSG00000157014", // TATDN2
//                "ENSG00000164362", // TERT
//                "ENSG00000163486", // SRGAP2
//                "ENSG00000242802" // AP5Z1
//        );
        copyLinesMatching(new File(fullEnsemblDataDir, "ensembl_gene_data.csv"), outputDir, geneIds);
        copyLinesMatching(new File(fullEnsemblDataDir, "ensembl_trans_amino_acids.csv"), outputDir, geneIds);
        copyLinesMatching(new File(fullEnsemblDataDir, "ensembl_trans_exon_data.csv"), outputDir, geneIds);
        copyLinesMatching(new File(fullEnsemblDataDir, "ensembl_trans_splice_data.csv"), outputDir, geneIds);
    }

    private static void copyLinesMatching(File input, File outputDir, Set<String> toFind) throws IOException
    {
        File outputFile = new File(outputDir, input.getName());
        outputFile.createNewFile();
        try (BufferedReader reader = new BufferedReader(new FileReader(input));
                BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {

            String line;
            boolean first = true; // Include the header
            while ((line = reader.readLine()) != null) {
                if (first) {
                    writer.write(line);
                    writer.newLine();
                    first = false;
                }
                if (matches(line, toFind)) {
                    writer.write(line);
                    writer.newLine();
                }
            }
        } catch (IOException e) {
            System.err.println("Error: " + e.getMessage());
        }
    }

    private static boolean matches(String line, Set<String> toFind)
    {
        if (toFind.isEmpty()) {
            return false;
        }
        for (String s : toFind)
        {
            if (line.contains(s))
            {
                return true;
            }
        }
        return false;
    }

//    @Test
    public void produceReducedChrFile() throws IOException
    {
        File outputDir = new File("/Users/timlavers/work/junk");

        RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")));
        int start = 4_000_000; //
        int end = start + 1_000_000;
        var chr = refGenomeSource.getBaseString("chr7", start, end);
//        System.out.println(chr.substring(10000, 10100));
        File chrFile = new File(outputDir, "chr7_part_4.txt");
        Files.writeString(chrFile.toPath(), chr, StandardCharsets.UTF_8);
    }
}
