package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.BamFragmentAllocator.calcFragmentLength;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.ReadCountsTest.REF_BASE_STR_1;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.isofox.TestUtils.TRANS_1;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createGeneReadData;
import static com.hartwig.hmftools.isofox.TestUtils.createReadRecord;
import static com.hartwig.hmftools.isofox.TestUtils.createRegion;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;
import static com.hartwig.hmftools.isofox.common.TransMatchType.ALT;
import static com.hartwig.hmftools.isofox.common.TransMatchType.EXONIC;
import static com.hartwig.hmftools.isofox.common.TransMatchType.SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.common.TransMatchType.UNSPLICED;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.common.FragmentType;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.results.ResultsWriter;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class TransClassificationTest
{
    @Test
    public void testReadTranscriptClassification()
    {
        int trans1 = 1;

        RegionReadData region = createRegion("GEN01", trans1, 1, "1", 100, 200);

        // unspliced read does not support the transcript
        ReadRecord read = createReadRecord(1, "1", 90, 110, REF_BASE_STR_1, createCigar(0, 21, 0));

        List<RegionReadData> regions = Lists.newArrayList(region);
        read.processOverlappingRegions(regions);

        assertEquals(UNSPLICED, read.getTranscriptClassification(trans1));

        // exonic read does support the transcript
        read = createReadRecord(1, "1", 120, 140, REF_BASE_STR_1, createCigar(0, 21, 0));

        regions = Lists.newArrayList(region);
        read.processOverlappingRegions(regions);

        assertEquals(EXONIC, read.getTranscriptClassification(trans1));

        // skipped exon
        RegionReadData region1 = createRegion("GEN01", trans1, 1, "1", 100, 120);
        RegionReadData region2 = createRegion("GEN01", trans1, 2, "1", 140, 160);
        RegionReadData region3 = createRegion("GEN01", trans1, 3, "1", 180, 200);

        // read covering multiple exons but skips the middle exon
        read = createReadRecord(1, "1", 110, 200, REF_BASE_STR_1, createCigar(0, 11, 59, 21, 0));

        regions = Lists.newArrayList(region1, region2, region3);
        read.processOverlappingRegions(regions);

        assertEquals(ALT, read.getTranscriptClassification(trans1));

        // incomplete intermediary exon
        Cigar cigar = new Cigar();
        cigar.add(new CigarElement(11, CigarOperator.M)); // matches half of first exon
        cigar.add(new CigarElement(19, CigarOperator.N));
        cigar.add(new CigarElement(11, CigarOperator.M)); // matches part of second exon
        cigar.add(new CigarElement(29, CigarOperator.N));
        cigar.add(new CigarElement(25, CigarOperator.M)); // matches past last exon into intron

        read = createReadRecord(1, "1", 110, 200, REF_BASE_STR_1, cigar);

        read.processOverlappingRegions(regions);

        assertEquals(ALT, read.getTranscriptClassification(trans1));

        // read spanning into intron (unspliced)
        cigar = new Cigar();
        cigar.add(new CigarElement(11, CigarOperator.M)); // matches half of first exon
        cigar.add(new CigarElement(19, CigarOperator.N));
        cigar.add(new CigarElement(31, CigarOperator.M)); // goes past second exon into intron
        cigar.add(new CigarElement(9, CigarOperator.N));
        cigar.add(new CigarElement(25, CigarOperator.M)); // matches past last exon into intron

        read = createReadRecord(1, "1", 110, 200, REF_BASE_STR_1, cigar);

        read.processOverlappingRegions(regions);

        assertEquals(UNSPLICED, read.getTranscriptClassification(trans1));

        // a read mapping to 2 regions, one as a splice junction, the other as exonic
        read = createReadRecord(1, "1", 181, 319, REF_BASE_STR_1 + REF_BASE_STR_1,
                createCigar(0, 20, 99, 20, 0));

        region1 = createRegion("GEN01", trans1, 1, "1", 100, 200);
        region2 = createRegion("GEN01", trans1, 2, "1", 300, 400);

        int trans2 = 2;
        region3 = createRegion("GEN01", trans2, 1, "1", 100, 220);

        regions = Lists.newArrayList(region1, region2, region3);
        read.processOverlappingRegions(regions);

        assertEquals(SPLICE_JUNCTION, read.getTranscriptClassification(trans1));
        assertEquals(ALT, read.getTranscriptClassification(trans2));
    }

    @Test
    public void testSoftClippedReadTranscriptClassification()
    {
        // test a read with soft-clipping matching the previous exon boundary
        RegionReadData region1 = createRegion(GENE_ID_1, TRANS_1, 1, "1", 100, 119);
        region1.setRefBases(REF_BASE_STR_1);
        RegionReadData region2 = createRegion(GENE_ID_1, TRANS_1, 2, "1", 200, 219);
        region2.setRefBases(REF_BASE_STR_1);
        RegionReadData region3 = createRegion(GENE_ID_1, TRANS_1, 3, "1", 300, 319);
        region3.setRefBases(REF_BASE_STR_1);

        region1.addPostRegion(region2);
        region2.addPreRegion(region1);
        region2.addPostRegion(region3);
        region3.addPreRegion(region2);

        String readBases = REF_BASE_STR_1.substring(REF_BASE_STR_1.length() - 1) + REF_BASE_STR_1 + REF_BASE_STR_1.substring(0, 10);
        ReadRecord read = createReadRecord(1, CHR_1, 200, 309, readBases, createCigar(1, 20, 80, 10, 0));

        List<RegionReadData> regions = Lists.newArrayList(region2, region3);
        read.processOverlappingRegions(regions);

        assertEquals(SPLICE_JUNCTION, read.getTranscriptClassification(TRANS_1));

        // any soft-clipped read which cannot be mapped is classified as ALT
        readBases = REF_BASE_STR_1 + "AAAAA";
        read = createReadRecord(1, CHR_1, 100, 119, readBases, createCigar(0, 20, 5));
        read.setFragmentInsertSize(200);
        read.processOverlappingRegions(Lists.newArrayList(region1));

        assertEquals(ALT, read.getTranscriptClassification(TRANS_1));

        readBases = "AAAAA" + REF_BASE_STR_1;
        read = createReadRecord(1, CHR_1, 300, 319, readBases, createCigar(5, 20, 0));
        read.setFragmentInsertSize(200);
        read.processOverlappingRegions(Lists.newArrayList(region3));

        assertEquals(ALT, read.getTranscriptClassification(TRANS_1));

        // likely adapter sequences are permitted
        read = createReadRecord(1, CHR_1, 300, 319, readBases, createCigar(5, 20, 0));
        read.setFragmentInsertSize(20);
        read.processOverlappingRegions(Lists.newArrayList(region3));

        assertEquals(EXONIC, read.getTranscriptClassification(TRANS_1));
    }

    @Test
    public void testFragmentLengthCalcs()
    {
        String transName1 = "TRANS01";

        TranscriptData transData = new TranscriptData(TRANS_1, transName1, GENE_NAME_1, true, (byte) 1,
                1000, 5000, null, null, "");

        transData.exons().add(new ExonData(TRANS_1, 1000, 1200, 1, -1, -1));
        transData.exons().add(new ExonData(TRANS_1, 2000, 2500, 2, -1, -1));
        transData.exons().add(new ExonData(TRANS_1, 4500, 5000, 3, -1, -1));
        transData.exons().add(new ExonData(TRANS_1, 5500, 6000, 4, -1, -1));

        // within 1st exon
        ReadRecord read1 = createReadRecord(1, CHR_1, 1010, 1029, REF_BASE_STR_1, createCigar(0, 20, 0));
        ReadRecord read2 = createReadRecord(1, CHR_1, 1170, 1199, REF_BASE_STR_1, createCigar(0, 20, 0));

        int fragLength = calcFragmentLength(transData, read1, read2);
        assertEquals(190, fragLength);

        // spanning 2 exons, both exonic
        read1 = createReadRecord(1, CHR_1, 1170, 1189, REF_BASE_STR_1, createCigar(0, 20, 0));
        read2 = createReadRecord(1, CHR_1, 2010, 2019, REF_BASE_STR_1, createCigar(0, 20, 0));

        fragLength = calcFragmentLength(transData, read1, read2);
        assertEquals(31 + 20, fragLength);

        // spanning 3 exons, both exonic
        read1 = createReadRecord(1, CHR_1, 1170, 1189, REF_BASE_STR_1, createCigar(0, 20, 0));
        read2 = createReadRecord(1, CHR_1, 4510, 4519, REF_BASE_STR_1, createCigar(0, 20, 0));

        fragLength = calcFragmentLength(transData, read1, read2);
        assertEquals(31 + 501 + 20, fragLength);

        // with 2 split reads
        read1 = createReadRecord(1, CHR_1, 1191, 2009, REF_BASE_STR_1, createCigar(0, 10, 799, 10, 0));
        read2 = createReadRecord(1, CHR_1, 2491, 4509, REF_BASE_STR_1, createCigar(0, 10, 1999, 10, 0));

        fragLength = calcFragmentLength(transData, read1, read2);
        assertEquals(10 + 501 + 10, fragLength);

        // with 2 split reads skipping an exon
        read1 = createReadRecord(1, CHR_1, 1191, 2009, REF_BASE_STR_1, createCigar(0, 10, 799, 10, 0));
        read2 = createReadRecord(1, CHR_1, 4991, 5509, REF_BASE_STR_1, createCigar(0, 10, 1999, 10, 0));

        fragLength = calcFragmentLength(transData, read1, read2);
        assertEquals(10 + 501 + 501 + 10, fragLength);
    }

    @Test
    public void testFragmentReadPairs()
    {
        // one read outside the gene
        IsofoxConfig config = new IsofoxConfig();
        config.Functions.clear();
        config.Functions.add(TRANSCRIPT_COUNTS);
        BamFragmentAllocator bamReader = new BamFragmentAllocator(config, new ResultsWriter(config));

        String transName1 = "TRANS01";

        TranscriptData transData1 = new TranscriptData(TRANS_1, transName1, GENE_NAME_1, true, (byte) 1,
                1000, 5000, null, null, "");

        transData1.exons().add(new ExonData(TRANS_1, 1000, 1200, 1, -1, -1));
        transData1.exons().add(new ExonData(TRANS_1, 2000, 2500, 2, -1, -1));
        transData1.exons().add(new ExonData(TRANS_1, 4500, 5000, 3, -1, -1));

        GeneReadData geneReadData = createGeneReadData(GENE_NAME_1, "1", (byte) 1, 1000, 5000);
        geneReadData.setTranscripts(Lists.newArrayList(transData1));

        ReadRecord read1 = createReadRecord(1, CHR_1, 100, 200, REF_BASE_STR_1, createCigar(0, 10, 0));
        ReadRecord read2 = createReadRecord(1, CHR_1, 1050, 1150, REF_BASE_STR_1, createCigar(0, 20, 0));

        List<ReadRecord> reads = Lists.newArrayList(read1, read2);

        GeneCollection geneSet = new GeneCollection(0, Lists.newArrayList(geneReadData));
        bamReader.processReadRecords(geneSet, reads);

        int[] geneCounts = geneSet.getCounts();
        assertEquals(1, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(CHIMERIC)]);

        // exon to intronic read
        read1 = createReadRecord(1, CHR_1, 1300, 1350, REF_BASE_STR_1, createCigar(0, 50, 0));
        read1.setFragmentInsertSize(1100);
        read2 = createReadRecord(1, CHR_1, 2300, 2400, REF_BASE_STR_1, createCigar(0, 100, 0));
        read2.setFragmentInsertSize(-1100);

        reads = Lists.newArrayList(read1, read2);
        bamReader.processReadRecords(geneSet, reads);

        assertEquals(2, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(FragmentType.UNSPLICED)]);

        // fully intronic
        read1 = createReadRecord(1, CHR_1, 2600, 2650, REF_BASE_STR_1, createCigar(0, 50, 0));
        read1.setFragmentInsertSize(300);
        read2 = createReadRecord(1, CHR_1, 2800, 2900, REF_BASE_STR_1, createCigar(0, 100, 0));
        read2.setFragmentInsertSize(-300);

        geneSet.clearCounts();

        reads = Lists.newArrayList(read1, read2);
        bamReader.processReadRecords(geneSet, reads);

        assertEquals(1, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(FragmentType.UNSPLICED)]);

        // alternative splicing - first from reads with splits
        geneSet.clearCounts();

        read1 = createReadRecord(1, CHR_1, 1050, 6100, REF_BASE_STR_1, createCigar(0, 50, 5000, 100, 0));
        read1.setFragmentInsertSize(500);
        read2 = createReadRecord(1, CHR_1, 3100, 3300, REF_BASE_STR_1, createCigar(0, 100, 0));
        read2.setFragmentInsertSize(-500);

        reads = Lists.newArrayList(read1, read2);
        bamReader.processReadRecords(geneSet, reads);

        assertEquals(1, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(CHIMERIC)]);

        int longInsertSize = config.MaxFragmentLength + 100;

        // long exon to exon read, treated as supporting
        geneSet.clearCounts();

        read1 = createReadRecord(1, CHR_1, 1050, 1099, REF_BASE_STR_1, createCigar(0, 50, 0));
        read1.setFragmentInsertSize(longInsertSize);
        read2 = createReadRecord(1, CHR_1, 2100, 2199, REF_BASE_STR_1, createCigar(0, 100, 0));
        read2.setFragmentInsertSize(-longInsertSize);

        reads = Lists.newArrayList(read1, read2);
        bamReader.processReadRecords(geneSet, reads);

        assertEquals(1, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(FragmentType.TRANS_SUPPORTING)]);

        // alt splicing - exon to exon read skipping an exon and long, currently not detected
        geneSet.clearCounts();

        read1 = createReadRecord(1, CHR_1, 1050, 1100, REF_BASE_STR_1, createCigar(0, 50, 0));
        read1.setFragmentInsertSize(longInsertSize);
        read2 = createReadRecord(1, CHR_1, 4550, 4650, REF_BASE_STR_1, createCigar(0, 100, 0));
        read2.setFragmentInsertSize(-longInsertSize);

        reads = Lists.newArrayList(read1, read2);
        bamReader.processReadRecords(geneSet, reads);

        assertEquals(1, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(FragmentType.ALT)]);
    }

}
