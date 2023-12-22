package com.hartwig.hmftools.svassembly;

import static org.assertj.core.api.Assertions.assertThat;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.Strings;
import com.hartwig.hmftools.svassembly.models.Record;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class ReadRescueTest
{
    private static final String REFERENCE_SEQUENCE = "/Users/james/code/data/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";

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
        assertThat(rescued).isNotSameAs(record);
        assertThat(rescued.getBasesString()).isEqualTo(record.getBasesString());
        assertThat(rescued.getBaseQuality()).isNotEqualTo(record.getBaseQuality());
        assertThat(rescued.getBaseQualityString()).isEqualTo("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF,F,:::FFFF,F:F:F,FF::FFFF,FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF,,,F,FFFFFFFF,F,F,FFFFFFFF,FF,,FF,FFFFFFF,FF,FFF");
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
        assertThat(rescued).isNotSameAs(record);
        assertThat(rescued.getBasesString()).isEqualTo(record.getBasesString());
        assertThat(rescued.getBaseQuality()).isNotEqualTo(record.getBaseQuality());
        assertThat(rescued.getBaseQualityString()).isEqualTo(Strings.reverseString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF,F,:::FFFF,F:F:F,FF::FFFF,FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"));
    }
}