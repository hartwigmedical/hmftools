package htsjdk.samtools

import junit.framework.Assert.assertEquals
import org.junit.Test

class SAMUtilsExtraTest {
    @Test
    fun testSoftClipsBeforeAndAfterMatch() {
        val cigarOperators = "10S".asCigar() + "130M".asCigar() + "10S".asCigar()
        val cigar = Cigar.fromCigarOperators(cigarOperators)!!

        val result = SAMUtilsExtra.getAlignmentBlocksWithSoftClips(cigar, 100, "TypeName")
        assertEquals(3, result.size)
        assertAlignment(result[0], 90, 1, 10)
        assertAlignment(result[1], 100, 11, 130)
        assertAlignment(result[2], 230, 141, 10)
    }

    @Test
    fun testSoftClipsInterruptedByIndelBeforeAndAfterMatch() {
        val cigarOperators = "10S".asCigar() + "1D".asCigar() + "130M".asCigar() + "1D".asCigar() + "10S".asCigar()
        val cigar = Cigar.fromCigarOperators(cigarOperators)!!

        val result = SAMUtilsExtra.getAlignmentBlocksWithSoftClips(cigar, 100, "TypeName")
        assertEquals(1, result.size)
        assertAlignment(result[0], 101, 11, 130)
    }

    private fun assertAlignment(victim: AlignmentBlock, expectedReferenceStart: Int, expectedReadStart: Int, expectedLength: Int) {
        assertEquals(expectedReferenceStart, victim.referenceStart)
        assertEquals(expectedReadStart, victim.readStart)
        assertEquals(expectedLength, victim.length)
    }

    private fun String.asCigar(): List<CigarOperator> {
        val op = CigarOperator.valueOf(this.takeLast(1))
        val count = this.take(this.length - 1).toInt()
        return cigar(count, op)
    }

    private fun cigar(count: Int, operator: CigarOperator): List<CigarOperator> {
        return (0 until count).map { operator }
    }

}