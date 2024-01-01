package com.hartwig.hmftools.esvee;

public class ReadRescueTest
{
    private static final String REFERENCE_SEQUENCE = "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";

    /* CHASHA FIXME
    @Test
    public void canRescueForwards() {
        final SAMRecord inner = new SAMRecord(new SAMFileHeader());
        final Record record = new Record(inner);
        inner.setReadString("TCACCCTCCTGAGTAGCTGGTGTGCACCATCACGCTCAGCTAATTTTTTTGGTTTTTTTTTGTTTTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTTTTGTTGGGGTTTAAACCTGTTGCCCCGGGGGGGCTCAAACACCGGGG");
        inner.setBaseQualityString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF,F,:::FFFF,F:F:F,FF::FFFF,FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF,F,,,,,:,::F,F,,:,,,,,,F:,,:,,:,,:F,,,,,::,,,,,,FF");
        inner.setReadName("A00624:8:HHKYHDSXX:3:1453:9968:28213");
        inner.setReferenceName("1");
        inner.setAlignmentStart(187626382);
        inner.setCigarString("61M90S");
        inner.setMappingQuality(60);

        final Record rescued = new ReadRescue(RefGenomeSource.loadRefGenome(REFERENCE_SEQUENCE)).rescueRead(record);
        assertTrue(rescued).isNotSameAs(record);
        assertTrue(rescued.getBasesString()).isEqualTo(record.getBasesString());
        assertTrue(rescued.getBaseQuality()).isNotEqualTo(record.getBaseQuality());
        assertTrue(rescued.getBaseQualityString()).isEqualTo("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF,F,:::FFFF,F:F:F,FF::FFFF,FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF,,,F,FFFFFFFF,F,F,FFFFFFFF,FF,,FF,FFFFFFF,FF,FFF");
    }

    @Test
    public void canRescueBackwards() {
        final SAMRecord inner = new SAMRecord(new SAMFileHeader());
        final Record record = new Record(inner);
        inner.setReadString("TCACCCTCCTGAGTAGCTGGTGTGCACCATCACGCTCAGCTAATTTTTTTGGTTTTTTTTTGTTTTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTTTTGTTGGGGTTTAAACCTGTTGCCCCGGGGGGGCTCAAACACCGGGG");
        inner.setBaseQualityString(Strings.reverseString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF,F,:::FFFF,F:F:F,FF::FFFF,FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF,F,,,,,:,::F,F,,:,,,,,,F:,,:,,:,,:F,,,,,::,,,,,,FF"));
        inner.setReadName("A00624:8:HHKYHDSXX:3:1453:9968:28213");
        inner.setReadNegativeStrandFlag(true);
        inner.setReferenceName("1");
        inner.setAlignmentStart(187626382 + 90);
        inner.setCigarString("90S61M");
        inner.setMappingQuality(60);

        final Record rescued = new ReadRescue(RefGenomeSource.loadRefGenome(REFERENCE_SEQUENCE)).rescueRead(record);
        assertTrue(rescued).isNotSameAs(record);
        assertTrue(rescued.getBasesString()).isEqualTo(record.getBasesString());
        assertTrue(rescued.getBaseQuality()).isNotEqualTo(record.getBaseQuality());
        assertTrue(rescued.getBaseQualityString()).isEqualTo(Strings.reverseString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF,F,:::FFFF,F:F:F,FF::FFFF,FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    }
    */
}