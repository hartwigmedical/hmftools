package com.hartwig.hmftools.pave.transval;

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

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class TestDataCreator
{

//    @Test
    public void createReducedEnsemblDataSet() throws IOException
    {
        File fullEnsemblDataDir = new File("/Users/timlavers/work/data/v6_0/ref/38/common/ensembl_data");
        File outputDir = new File("/Users/timlavers/work/junk/ensemblmini");

        final Set<String> geneNames = Set.of(
                "MTOR,",
                "TET2,",
                "VHL,",
                "BRAF,",
                "ADCK2,",
                "EGFR,",
                "RNU1-82P,",
                "ZYX,",
                "PIK3R1"
        );
        copyLinesMatching(new File(fullEnsemblDataDir, "ensembl_gene_data.csv"), outputDir, geneNames);
        copyLinesMatching(new File(fullEnsemblDataDir, "ensembl_protein_features.csv"), outputDir, Set.of());
        Set<String> geneIds = Set.of(
                "ENSG00000198793", // MTOR
                "ENSG00000168769", // TET2
                "ENSG00000134086",  // VHL
                "ENSG00000133597", // ADCK2
                "ENSG00000157764", // BRAF
                "ENSG00000146648", // EGFR
                "ENSG00000212153", // RNU-82P
                "ENSG00000159840", // ZYX
                "ENSG00000145675" // PIK3R1
        );
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
        var chromosomeLengths = refGenomeSource.chromosomeLengths();
        System.out.println(chromosomeLengths.size());
        int chrLength = chromosomeLengths.get("chr5");
        System.out.println(chrLength);
        int start = 68_000_000; //
        int end = start + 1_000_000;
        var chr = refGenomeSource.getBaseString("chr5", start, end);
//        System.out.println(chr.substring(10000, 10100));
        File chrFile = new File(outputDir, "chr5_part_68.txt");
        Files.writeString(chrFile.toPath(), chr, StandardCharsets.UTF_8);
    }

//    @Test
    public void braf() throws IOException
    {
        File ensemblDataDir = new File("/Users/timlavers/work/data/v6_0/ref/38/common/ensembl_data");
        File ensembleGeneDataFile = new File(ensemblDataDir, "ensembl_gene_data.csv");

        EnsemblDataCache ensembl = new EnsemblDataCache(ensemblDataDir.getAbsolutePath(), RefGenomeVersion.V38);
        ensembl.setRequiredData(true, true, true, true);
        ensembl.load(false);
        ensembl.createTranscriptIdMap();
        GeneData brafData = ensembl.getGeneDataByName("BRAF");
        Assert.assertEquals("chr7", brafData.Chromosome);

        List<TranscriptData> brafTranscriptData = ensembl.getTranscripts(brafData.GeneId);
        Assert.assertEquals(1, brafTranscriptData.size()); // Is it just always loading the canonical one?

        TranscriptData canonicalBraf = brafTranscriptData.get(0);
        Assert.assertEquals(18, canonicalBraf.exons().size());
        ExonData exon0 = canonicalBraf.exons().get(0);
//        Assert.assertEquals(140734597, exon0.Start);
        Assert.assertEquals(140734770, exon0.End);

        ExonData exon1 = canonicalBraf.exons().get(1);
        Assert.assertEquals(140739812, exon1.Start);
        Assert.assertEquals(140739946, exon1.End);


        RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")));
//        String bases0 = refGenomeSource.getBaseString(brafData.Chromosome, exon0.Start, exon0.End);
//        System.out.println(bases0);
        String bases1 = refGenomeSource.getBaseString(brafData.Chromosome, exon1.Start, exon1.End);
        String basesFromCCDS = "ATAATTTTTATGGTGGGACGAGGATACCTGTCTCCAGA"
                + "TCTCAGTAAGGTACGGAGTAACTGTCCAAAAGCCATGAAGAGATTAATGGCAGAGTGCCTCAAAAAGAAA"
                + "AGAGATGAGAGACCACTCTTTCCCCAA";
        Assert.assertEquals(basesFromCCDS, Nucleotides.reverseComplementBases(bases1));
//        System.out.println(bases1);

        // Data for BRAF p.Val600Glu (p.V600E)
        ExonData exon3 = canonicalBraf.exons().get(3);
        Assert.assertEquals(140753275, exon3.Start);
        Assert.assertEquals(140753393, exon3.End);

        String bases3 = refGenomeSource.getBaseString(brafData.Chromosome, exon3.Start, exon3.End);
        String bases3FromCCDS = "ATATATTTCTTCATGAAGACCTCACAGTAAAAATAGGTGATTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGATCCATTTTGTGGATG";
        Assert.assertEquals(bases3FromCCDS, Nucleotides.reverseComplementBases(bases3));

        TinyGenome tinyGenome = new TinyGenome();
        String bases3FromChr7Genome = tinyGenome.getBaseString(brafData.Chromosome, exon3.Start, exon3.End);
        Assert.assertEquals(bases3, bases3FromChr7Genome);

        // Codon for V at 600 is GTG
        // 600-602 is VKS = GTGAAATCT
        // For V600E need GTG -> GAA or GTG -> GAG. The single nucleotide change GTG -> GAG is most likely.
        // The position in the coding sequence of this is 1799 (-1 as it's on the reverse strand).
        // So the output is c.

//        Map<String, TranscriptAminoAcids> transAminoAcidMap = Maps.newHashMap();
//        EnsemblDataLoader.loadTranscriptAminoAcidData(
//                ensemblDataDir, transAminoAcidMap, Lists.newArrayList(), false);
//
//        var mapSize = transAminoAcidMap.size();
//        System.out.println(mapSize);
//
//        var chromosomeLengths = refGenomeSource.chromosomeLengths();
//        System.out.println(chromosomeLengths.size());
//        int chr7Length = chromosomeLengths.get("chr7");
//        int start = 140_000_000;
//        int end = start + 10_000_000;
//        var chr7 = refGenomeSource.getBaseString("chr7", start, end);
//        System.out.println(chr7.substring(10000, 10100));
//        File chr7File = new File(ensemblDataDir, "chr7_part.txt");
//        Files.writeString(chr7File.toPath(), chr7, StandardCharsets.UTF_8);
    }
}
