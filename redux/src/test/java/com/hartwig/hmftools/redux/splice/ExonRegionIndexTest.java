package com.hartwig.hmftools.redux.splice;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.Test;

// contains() uses binary search on merged ranges; earlier scan-back was O(N) and tripped on overlap edge cases.
public class ExonRegionIndexTest
{
    @Test
    public void containsHitsIncludeBoundaries() throws Exception
    {
        ExonRegionIndex idx = build(List.of(new int[] { 100, 200 }, new int[] { 300, 400 }));

        assertTrue(idx.contains("chr1", 100));
        assertTrue(idx.contains("chr1", 150));
        assertTrue(idx.contains("chr1", 200));
        assertTrue(idx.contains("chr1", 300));
        assertTrue(idx.contains("chr1", 400));
    }

    @Test
    public void containsMissesIntronicAndIntergenic() throws Exception
    {
        ExonRegionIndex idx = build(List.of(new int[] { 100, 200 }, new int[] { 300, 400 }));

        assertFalse(idx.contains("chr1", 50));
        assertFalse(idx.contains("chr1", 250));
        assertFalse(idx.contains("chr1", 500));
    }

    @Test
    public void containsAcrossOverlappingAndAbuttingIntervals() throws Exception
    {
        ExonRegionIndex idx = build(List.of(
                new int[] { 100, 200 },
                new int[] { 150, 300 },
                new int[] { 301, 400 },
                new int[] { 600, 700 }));

        assertTrue(idx.contains("chr1", 250));
        assertTrue(idx.contains("chr1", 301));
        assertTrue(idx.contains("chr1", 400));
        assertFalse(idx.contains("chr1", 401));
        assertFalse(idx.contains("chr1", 599));
        assertTrue(idx.contains("chr1", 700));
    }

    @Test
    public void chromosomeNormalizationIsBidirectional() throws Exception
    {
        ExonRegionIndex idx = build(List.of(new int[] { 100, 200 }));
        assertTrue(idx.contains("chr1", 150));
        assertTrue(idx.contains("1", 150));
        assertFalse(idx.contains("chr2", 150));
    }

    private static ExonRegionIndex build(final List<int[]> exons) throws Exception
    {
        final Path dir = Files.createTempDirectory("exonIdxTest");
        Files.writeString(dir.resolve("ensembl_gene_data.csv"),
                "GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd\nENSG_TEST,TESTG,1,1,1,100000\n");
        final StringBuilder sb = new StringBuilder("GeneId,CanonicalTranscriptId,Strand,TransId,TransName,BioType,"
                + "TransStart,TransEnd,ExonRank,ExonStart,ExonEnd,ExonPhase,ExonEndPhase,CodingStart,CodingEnd,RefSeqId\n");
        int rank = 1;
        for(int[] ex : exons)
        {
            sb.append("ENSG_TEST,1,1,1,ENST_TEST,protein_coding,1,100000,").append(rank++).append(",")
                    .append(ex[0]).append(",").append(ex[1]).append(",0,0,0,0,\n");
        }
        Files.writeString(dir.resolve("ensembl_trans_exon_data.csv"), sb.toString());
        return ExonRegionIndex.load(dir.toString());
    }
}
